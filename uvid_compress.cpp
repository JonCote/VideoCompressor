/* uvid_compress.cpp

    Jonathan Cote, V00962634
*/
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <utility>
#include <queue>
#include <map>
#include <bitset>
#include <unordered_map>
#include <tuple>
#include <math.h>
#include <limits>
#include "output_stream.hpp"
#include "yuv_stream.hpp"
#include "uvid_common.hpp"



std::vector<std::vector<int>> Lumi = 
    {
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55},
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };
std::vector<std::vector<int>> Chromi = 
    {
        {17, 18, 24, 47, 99, 99, 99, 99},
        {18, 21, 26, 66, 99, 99, 99, 99},
        {24, 26, 56, 99, 99, 99, 99, 99},
        {47, 66, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99}
    };


std::vector<std::vector<double>> C;
std::vector<std::vector<double>> C_T;


void make_C() {
    C = create_2d_vector<double>(8, 8);
    C_T = create_2d_vector<double>(8, 8);
    for(unsigned int i = 0; i < 8; ++i) {
        for(unsigned int j = 0; j < 8; ++j) {
            double C_i;
            if(i == 0){
                C_i = 1.0 / sqrt(8.0);
            }else {
                C_i = sqrt(2.0 / 8.0);
            }
            C.at(i).at(j) = C_i * (cos(((2.0 * j + 1.0) * i * M_PI) / 16.0));
            
            C_T.at(j).at(i) = C.at(i).at(j);
        }
    }
}


//block_DCT(Y_block, Y_curr, Y_prev, frame_out, local_motion_vector.at(0), local_motion_vector.at(1), height, width, row, col, false);
void block_DCT(std::vector<std::vector<int>>& block, std::vector<std::vector<int>>& decomp,
                std::vector<std::vector<int>>& prev_frame, std::vector<int>& out_stream,
                unsigned int prev_start_x, unsigned int prev_start_y, unsigned int height,
                unsigned int width, int block_start_y, int block_start_x, bool use_chromi,
                unsigned int frame_type)
{

    // subdivide block into 8x8 sections (input block must be a NxN with N % 8 = 0)
    unsigned int block_size = block.size();
    auto temp_prev = create_2d_vector<int>(block_size, block_size);
    for(unsigned int r = 0; r < block_size; r+=8) {
        for(unsigned int c = 0; c < block_size; c+=8) {
            
            if(use_chromi) {
                calc_Block_Quantized_DCT(block, block, Chromi, C, C_T, c, r);
            }
            else {
                calc_Block_Quantized_DCT(block, block, Lumi, C, C_T, c, r);
            }

            // decomp for P-Frame
            if(use_chromi) {
                revs_Block_Quantized_DCT(block, temp_prev, Chromi, C, C_T, c, r);
            }
            else {
                revs_Block_Quantized_DCT(block, temp_prev, Lumi, C, C_T, c, r);
            }

            for(unsigned int i = 0; i < 64; ++i) {
                out_stream.push_back(block.at(Zigzag.at(i).first + r).at(Zigzag.at(i).second + c));
            }
            
        }
    }

    for(unsigned int i = 0; i < block_size; ++i){
        for(unsigned int j = 0; j < block_size; ++j) {
            if(block_start_y + i < height && block_start_x + j < width) {
                if(frame_type == 1) {
                    unsigned int prev_frame_row = prev_start_y + i;
                    unsigned int prev_frame_col = prev_start_x + j;

                    if(prev_frame_row >= height) {
                        prev_frame_row = height - 1;  // Copy last row
                    }
                    if(prev_frame_col >= width) {
                        prev_frame_col = width - 1;  // Copy last col
                    }

                    decomp.at(block_start_y + i).at(block_start_x + j) = temp_prev.at(i).at(j) + prev_frame.at(prev_frame_row).at(prev_frame_col);
                }
                else {
                    decomp.at(block_start_y + i).at(block_start_x + j) = temp_prev.at(i).at(j);
                }
                
            }
        } 
    }
}


float calc_AAD(std::vector<std::vector<int>>& block, std::vector<std::vector<int>>& prev_frame, unsigned int height, unsigned int width, unsigned int dx, unsigned int dy) {
    float aad = 0.0;

    for(int y = 0; y < 16; ++y) {
        for(int x = 0; x < 16; ++x) {
            unsigned int frame_row = y + dy;
            unsigned int frame_col = x + dx;

            if(frame_row >= height) {
                frame_row = height - 1;  // Copy last row
            }
            if(frame_col >= width) {
                frame_col = width - 1;  // Copy last col
            }

            aad += std::abs(block.at(y).at(x) - prev_frame.at(frame_row).at(frame_col));
        }
    }

    return aad / 256.0;
}



std::vector<unsigned int> binary_motion_vec_search(std::vector<std::vector<int>>& block, std::vector<std::vector<int>>& prev_frame,
                                                  unsigned int block_y_start, unsigned int block_x_start, unsigned int height, unsigned int width)
{
    // binary search
    unsigned int dx = block_x_start;
    unsigned int dy = block_y_start;
    float min_AAD = calc_AAD(block, prev_frame, height, width, dx, dy);
    if(min_AAD < 5.0) {
        return std::vector<unsigned int>{dx, dy};
    }

    unsigned int curr_y_start = 0;
    unsigned int curr_x_start = 0;
    unsigned int curr_y_end = height;
    unsigned int curr_x_end = width;
    unsigned int curr_dx = width/2 - 8;
    unsigned int curr_dy = height/2 - 8;

    float aad = calc_AAD(block, prev_frame, height, width, curr_dx, curr_dy);
    if(aad < min_AAD) {
        min_AAD = aad;
        dx = curr_dx;
        dy = curr_dy;
    }

    while(curr_y_end - curr_y_start >= 32 && curr_x_end - curr_x_start >= 32) {
        float temp_aad = std::numeric_limits<float>::max();
        float local_low_aad = std::numeric_limits<float>::max();
        
        // spread out evenly in 4 directions
        // get start and end of sub boxes in 4 directions
        unsigned int y_start_1 = curr_y_start;
        unsigned int x_start_1 = curr_x_start;
        unsigned int y_end_1 = curr_y_end - ((curr_y_end - curr_y_start)/2);
        unsigned int x_end_1 = curr_x_end - ((curr_x_end - curr_x_start)/2);

        unsigned int y_start_2 = curr_y_start;
        unsigned int x_start_2 = curr_x_start + ((curr_x_end - curr_x_start)/2);
        unsigned int y_end_2 = curr_y_end - ((curr_y_end - curr_y_start)/2);
        unsigned int x_end_2 = curr_x_end;

        unsigned int y_start_3 = curr_y_start + ((curr_y_end - curr_y_start)/2);
        unsigned int x_start_3 = curr_x_start;
        unsigned int y_end_3 = curr_y_end;
        unsigned int x_end_3 = curr_x_end - ((curr_x_end - curr_x_start)/2);

        unsigned int y_start_4 = curr_y_start + ((curr_y_end - curr_y_start)/2);
        unsigned int x_start_4 = curr_x_start + ((curr_x_end - curr_x_start)/2);
        unsigned int y_end_4 = curr_y_end;
        unsigned int x_end_4 = curr_x_end;

        // test center of each sub box and assign start, end values for best to curr
        unsigned int dy_1 = y_start_1 + (y_end_1 - y_start_1)/2 - 8;
        unsigned int dx_1 = x_start_1 + (x_end_1 - x_start_1)/2 - 8;

        unsigned int dy_2 = y_start_2 + (y_end_2 - y_start_2)/2 - 8;
        unsigned int dx_2 = x_start_2 + (x_end_2 - x_start_2)/2 - 8;

        unsigned int dy_3 = y_start_3 + (y_end_3 - y_start_3)/2 - 8;
        unsigned int dx_3 = x_start_3 + (x_end_3 - x_start_3)/2 - 8;

        unsigned int dy_4 = y_start_4 + (y_end_4 - y_start_4)/2 - 8;
        unsigned int dx_4 = x_start_4 + (x_end_4 - x_start_4)/2 - 8;


        temp_aad = calc_AAD(block, prev_frame, height, width, dx_1, dy_1);
        if(temp_aad <= local_low_aad) {
            local_low_aad = temp_aad;
            curr_y_start = y_start_1;
            curr_x_start = x_start_1;
            curr_y_end = y_end_1;
            curr_x_end = x_end_1;
            curr_dx = dx_1;
            curr_dy = dy_1;
        }

        temp_aad = calc_AAD(block, prev_frame, height, width, dx_2, dy_2);
        if(temp_aad <= local_low_aad) {
            local_low_aad = temp_aad;
            curr_y_start = y_start_2;
            curr_x_start = x_start_2;
            curr_y_end = y_end_2;
            curr_x_end = x_end_2;
            curr_dx = dx_2;
            curr_dy = dy_2;
        }

        temp_aad = calc_AAD(block, prev_frame, height, width, dx_3, dy_3);
        if(temp_aad <= local_low_aad) {
            local_low_aad = temp_aad;
            curr_y_start = y_start_3;
            curr_x_start = x_start_3;
            curr_y_end = y_end_3;
            curr_x_end = x_end_3;
            curr_dx = dx_3;
            curr_dy = dy_3;
        }

        temp_aad = calc_AAD(block, prev_frame, height, width, dx_4, dy_4);
        if(temp_aad <= local_low_aad) {
            local_low_aad = temp_aad;
            curr_y_start = y_start_4;
            curr_x_start = x_start_4;
            curr_y_end = y_end_4;
            curr_x_end = x_end_4;
            curr_dx = dx_4;
            curr_dy = dy_4;
        }

        // assign dx, dy if aad < min_AAD
        if(local_low_aad < min_AAD) {
            min_AAD = local_low_aad;
            dx = curr_dx;
            dy = curr_dy;
        }

    }

    return std::vector<unsigned int>{dx, dy};
}


// adapted from provided arithmatic encoder by B. Bird
void arith_compress(std::vector<u32>& stream, std::vector<int>& data_stream) {
    const u32 EOF_SYMBOL = 255; //1027; //4095;

    std::array<u32, EOF_SYMBOL+1> frequencies {};
    frequencies.fill(1);

    std::array<u64, EOF_SYMBOL+2> CF_low {};
    CF_low.at(0) = 0;
    for (unsigned int i = 1; i < EOF_SYMBOL+2; i++){
        CF_low.at(i) = CF_low.at(i-1) + frequencies.at(i-1);
    }

    u64 global_cumulative_frequency = CF_low.at(EOF_SYMBOL+1);

    assert(global_cumulative_frequency <= 0xffffffff); //If this fails, frequencies must be scaled down


    u32 lower_bound = 0;  //Bit sequence of all zeros
    u32 upper_bound = ~0; //Bit sequence of all ones

    int underflow_counter = 0;

    for(std::size_t ind = 0; ind < data_stream.size(); ++ind){
        u32 symbol = data_stream.at(ind) + 127; //513; //+ 2047
        
        //For safety, we will use u64 for all of our intermediate calculations.
        u64 current_range = ((u64)upper_bound + 1) - (u64)lower_bound;
        u64 symbol_range_low = CF_low.at(symbol);
        u64 symbol_range_high = CF_low.at(symbol+1);
        upper_bound = lower_bound + (current_range*symbol_range_high)/global_cumulative_frequency - 1;
        lower_bound = lower_bound + (current_range*symbol_range_low)/global_cumulative_frequency;

        // <-- This is probably where we would adjust the frequency table if we used an adaptive model.
        // update frequency table
        frequencies.at(symbol)++;
        
        for(unsigned int i = symbol + 1; i < EOF_SYMBOL+2; i++) {
            CF_low.at(i) = CF_low.at(i-1) + frequencies.at(i-1);
        }

        global_cumulative_frequency = CF_low.at(EOF_SYMBOL+1);
        assert(global_cumulative_frequency <= 0xffffffff);


        //Now determine if lower_bound and upper_bound share any of their most significant bits and push
        //them to the output stream if so.

        while(1){
            //Check if most significant bits (bit index 31) match.
            if ((upper_bound>>31) == (lower_bound>>31)){ 
                //Push the most significant bit of upper/lower
                u32 b = (upper_bound>>31);
                stream.push_back(b);
                //Now push underflow_counter copies of the opposite bit
                for(int i = 0; i < underflow_counter; i++){
                    stream.push_back(!b);
                }
                underflow_counter = 0;

                //Shift out the MSB of upper_bound (and shift in a 1 from the right)
                upper_bound <<= 1;
                upper_bound |= 1;

                //Shift out the MSB of lower_bound (and allow a 0 to be shifted in from the right)
                lower_bound <<= 1;
                
            }else if ( ((lower_bound>>30)&0x1) == 1 && ((upper_bound>>30)&0x1) == 0){
                //If the MSBs didn't match, then the MSB of upper_bound must be 1 and
                //the MSB of lower_bound must be 0.
                //If we discover that lower_bound = 01... and upper_bound = 10... 
                //(which is what the if-statement above tests), then we have
                //to account for underflow.

                underflow_counter++;

                //If upper_bound = 10(xyz...), set upper_bound = 1(xyz...)
                //(that is, splice out the second-most-significant bit)
                upper_bound <<= 1;
                upper_bound |= (1U<<31);
                upper_bound |= 1;

                //If lower_bound = 01(abc...), set lower_bound = 0(abd...)
                lower_bound <<= 1;
                lower_bound &= (1U<<31) - 1; //i.e. 0x7fffffff

            }else{
                break;
            }
        }
    }
    stream.push_back(0);
    stream.push_back(1);

    // Emit 1s to fill out the byte
    std::size_t byte_off = 8 - stream.size() % 8;
    for(std::size_t i = 0; i<byte_off; ++i) {
        stream.push_back(1);
    }
}


// lossless compression stuff end


int main(int argc, char** argv){

    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }
    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};

    // quality setting
    unsigned int quality_code = 2;  // default medium quality
    if(quality.compare("low") == 0) {
        for(unsigned int i = 0; i < 8; ++i){
            for(unsigned int j = 0; j < 8; ++j) {
                Lumi.at(i).at(j) = Lumi.at(i).at(j) * 2;
                Chromi.at(i).at(j) = Chromi.at(i).at(j) * 2;
            }
        }
        quality_code = 1;
    }else if(quality.compare("high") == 0) {
        for(unsigned int i = 0; i < 8; ++i){
            for(unsigned int j = 0; j < 8; ++j) {
                Lumi.at(i).at(j) = Lumi.at(i).at(j) / 2;
                Chromi.at(i).at(j) = Chromi.at(i).at(j) / 2;
            }
        }
        quality_code = 3;
    }
    else {
        for(unsigned int i = 0; i < 8; ++i){
            for(unsigned int j = 0; j < 8; ++j) {
                Lumi.at(i).at(j) = Lumi.at(i).at(j);
                Chromi.at(i).at(j) = Chromi.at(i).at(j);
            }
        }
        quality_code = 2;
    }
    make_C();

    YUVStreamReader reader {std::cin, width, height};
    OutputBitStream output_stream {std::cout};

    output_stream.push_u16(quality_code);
    output_stream.push_u32(height);
    output_stream.push_u32(width);

    auto Y_prev = create_2d_vector<int>(height,width);
    auto Cb_prev = create_2d_vector<int>(height/2,width/2);
    auto Cr_prev = create_2d_vector<int>(height/2,width/2);

    int frame_count = 0;

    while (reader.read_next_frame()){
        // read next frame
        output_stream.push_byte(1);  //Use a one byte flag to indicate whether there is a frame here

        YUVFrame420& frame = reader.frame();

        std::vector<int> frame_out {};
        unsigned int frame_type = 0;

        if(frame_count == 0) {
            output_stream.push_bit(0);  // one byte flag to indicate I-Frame
        }
        else {
            output_stream.push_bit(1); // one Byte flag to indicate P-Frame
            frame_type = 1;
        }

        auto Y_curr = create_2d_vector<int>(height,width);
        auto Cb_curr = create_2d_vector<int>(height/2,width/2);
        auto Cr_curr = create_2d_vector<int>(height/2,width/2);
        std::vector<std::vector<unsigned int>> motion_vectors {};

        for(unsigned int row = 0; row < height; row+=16) {
            for (unsigned int col = 0; col < width; col+=16) {
                //std::pair<int, int> block_inds;
                auto Y_block = create_2d_vector<int>(16, 16);
                std::vector<std::vector<int>> Cb_block(8, std::vector<int>(8, 0));
                std::vector<std::vector<int>> Cr_block(8, std::vector<int>(8, 0));

                // read a yuv 4:2:0 macroblocks
            
                // 16x16 block for Y
                for(unsigned int y = 0; y < 16; ++y) {
                    for(unsigned int x = 0; x < 16; ++x) {
                        unsigned int frame_row = row + y;
                        unsigned int frame_col = col + x;

                        if(frame_row >= height) {
                            frame_row = height - 1;  // Copy last row
                        }
                        if(frame_col >= width) {
                            frame_col = width - 1;  // Copy last col
                        }

                        Y_block.at(y).at(x) = frame.Y(frame_col, frame_row) - 128;
                    }
                }

                // make 1 8x8 block for Cb and Cr
                for(unsigned int y = 0; y < 8; ++y) {
                    for(unsigned int x = 0; x < 8; ++x) {
                        unsigned int frame_row = row/2 + y;
                        unsigned int frame_col = col/2 + x;

                        if(frame_row >= height/2) {
                            frame_row = height/2 - 1;  // Copy last row
                        }
                        if(frame_col >= width/2) {
                            frame_col = width/2 - 1;  // Copy last col
                        }

                        Cb_block.at(y).at(x) = frame.Cb(frame_col, frame_row) - 128;
                        Cr_block.at(y).at(x) = frame.Cr(frame_col, frame_row) - 128;
                    }
                }

                std::vector<unsigned int> local_motion_vector {col, row};

                if(frame_count != 0) {
                    // P-Frame

                    local_motion_vector = binary_motion_vec_search(Y_block, Y_prev, row, col, height, width);
                    
                    // subtract prev blocks for curr blocks
                    for(unsigned int y = 0; y < 8; ++y) {
                        for(unsigned int x = 0; x < 8; ++x) {
                            unsigned int frame_row = local_motion_vector.at(1)/2 + y;
                            unsigned int frame_col = local_motion_vector.at(0)/2 + x;

                            if(frame_row >= height/2) {
                                frame_row = height/2 - 1;  // Copy last row
                            }
                            if(frame_col >= width/2) {
                                frame_col = width/2 - 1;  // Copy last col
                            }

                            Cb_block.at(y).at(x) = Cb_block.at(y).at(x) - Cb_prev.at(frame_row).at(frame_col);
                            Cr_block.at(y).at(x) = Cr_block.at(y).at(x) - Cr_prev.at(frame_row).at(frame_col);
                        }
                    }
                    
                    for(int y = 0; y < 16; ++y) {
                        for(int x = 0; x < 16; ++x) {
                            unsigned int frame_row = local_motion_vector.at(1) + y;
                            unsigned int frame_col = local_motion_vector.at(0) + x;

                            if(frame_row >= height) {
                                frame_row = height - 1;  // Copy last row
                            }
                            if(frame_col >= width) {
                                frame_col = width - 1;  // Copy last col
                            }

                            Y_block.at(y).at(x) = Y_block.at(y).at(x) - Y_prev.at(frame_row).at(frame_col);
                        }
                    }
                    
                    // push motion vector to frame output
                    motion_vectors.push_back(local_motion_vector);
                }

                block_DCT(Y_block, Y_curr, Y_prev, frame_out, local_motion_vector.at(0), local_motion_vector.at(1), height, width, row, col, false, frame_type);
                block_DCT(Cb_block, Cb_curr, Cb_prev, frame_out, local_motion_vector.at(0)/2, local_motion_vector.at(1)/2, height/2, width/2, row/2, col/2, true, frame_type);
                block_DCT(Cr_block, Cr_curr, Cr_prev, frame_out, local_motion_vector.at(0)/2, local_motion_vector.at(1)/2, height/2, width/2, row/2, col/2, true, frame_type);

            }
        }

        if(frame_count == 10) {
            frame_count = 0;
        }
        else {
            frame_count++;
        }

        Y_prev = Y_curr;
        Cb_prev = Cb_curr;
        Cr_prev = Cr_curr;

        // push number of motion vectors for frame
        output_stream.push_u32(motion_vectors.size());

        // push the motion vectors
        for(unsigned int i = 0; i < motion_vectors.size(); ++i) {
            output_stream.push_u32(motion_vectors.at(i).at(0));
            output_stream.push_u32(motion_vectors.at(i).at(1));
        }

        // bitstream compression
        frame_out = rle_encode(frame_out);
        
        std::vector<u32> encoded_stream {};
        arith_compress(encoded_stream, frame_out);
        
        output_stream.push_u32(frame_out.size());
        output_stream.push_u32(encoded_stream.size());
        for(std::size_t i = 0; i < encoded_stream.size(); ++i) {
            output_stream.push_bit(encoded_stream.at(i));
        }
        
    }

    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();

    return 0;
}