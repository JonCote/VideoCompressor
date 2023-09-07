/* uvid_decompress.cpp
   
   Jonathan Cote, V00962634

*/
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include "input_stream.hpp"
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


void rev_block_DCT(std::vector<std::vector<int>>& block, bool use_chromi)
{
    unsigned int block_size = block.size();
    for(unsigned int r = 0; r < block_size; r+=8) {
        for(unsigned int c = 0; c < block_size; c+=8) {
            if(use_chromi) {
                revs_Block_Quantized_DCT(block, block, Chromi, C, C_T, c, r);
            }
            else {
                revs_Block_Quantized_DCT(block, block, Lumi, C, C_T, c, r);
            }
        }
    }
}



void place_P_frame_block(std::vector<std::vector<int>>& block, std::vector<std::vector<int>>& prev_frame, std::vector<std::vector<int>>& curr_frame, YUVFrame420& frame,
                            unsigned int dx, unsigned int dy, unsigned int height, unsigned int width,
                            unsigned int block_start_y, unsigned int block_start_x, unsigned int channel) 
{
    unsigned int block_size = block.size();
    for(unsigned int i = 0; i < block_size; ++i){
        for(unsigned int j = 0; j < block_size; ++j) {
            if(block_start_y + i < height && block_start_x + j < width) {
                unsigned int prev_frame_row = dy + i;
                unsigned int prev_frame_col = dx + j;

                if(prev_frame_row >= height) {
                    prev_frame_row = height - 1;  // Copy last row
                }
                if(prev_frame_col >= width) {
                    prev_frame_col = width - 1;  // Copy last col
                }

                int value = block.at(i).at(j) + prev_frame.at(prev_frame_row).at(prev_frame_col);

                // output to frame
                if(channel == 0) {
                    frame.Y(block_start_x + j, block_start_y + i) = round_and_clamp_to_char(value + 128);
                }
                if(channel == 1) {
                    frame.Cb(block_start_x + j, block_start_y + i) = round_and_clamp_to_char(value + 128);
                }
                if(channel == 2) {
                    frame.Cr(block_start_x + j, block_start_y + i) = round_and_clamp_to_char(value + 128);
                }

                // store for previous frame calculations
                curr_frame.at(block_start_y + i).at(block_start_x + j) = value;

            }
        } 
    }  
}

void place_I_frame_block(std::vector<std::vector<int>>& block, std::vector<std::vector<int>>& curr_frame, YUVFrame420& frame,
                         unsigned int height, unsigned int width, unsigned int block_start_y, unsigned int block_start_x, unsigned int channel) 
{
    unsigned int block_size = block.size();
    for(unsigned int i = 0; i < block_size; ++i){
        for(unsigned int j = 0; j < block_size; ++j) {
            if(block_start_y + i < height && block_start_x + j < width) {
                int value = block.at(i).at(j);
                
                // output to frame
                if(channel == 0) {
                    frame.Y(block_start_x + j, block_start_y + i) = round_and_clamp_to_char(value + 128);
                }
                if(channel == 1) {
                    frame.Cb(block_start_x + j, block_start_y + i) = round_and_clamp_to_char(value + 128);
                }
                if(channel == 2) {
                    frame.Cr(block_start_x + j, block_start_y + i) = round_and_clamp_to_char(value + 128);
                }

                // store for previous frame calculations
                curr_frame.at(block_start_y + i).at(block_start_x + j) = value;
                
            }
        } 
    } 
}


// lossless compression stuff starts

// adapted from provided arithmatic decoder by B. Bird
void arith_decompress(InputBitStream& stream, std::vector<int>& frame_data, unsigned int frame_len, unsigned int encoded_len) {
    const u32 EOF_SYMBOL = 255; //4095;

    // Create a static frequency table with a frequency of 1 for 
    // all symbols except lowercase/uppercase letters (symbols 65-122)

    std::array<u32, EOF_SYMBOL+1> frequencies {};
    frequencies.fill(1);

    //Now compute cumulative frequencies for each symbol.
    //We actually want the range [CF_low,CF_high] for each symbol,
    //but since CF_low(i) = CF_high(i-1), we only really have to compute
    //the array of lower bounds.

    //The cumulative frequency range for each symbol i will be 
    //[ CF_low.at(i), CF_low.at(i+1) ) 
    //(note that it's a half-open interval)
    std::array<u64, EOF_SYMBOL+2> CF_low {};
    CF_low.at(0) = 0;
    for (unsigned int i = 1; i < EOF_SYMBOL+2; i++){
        CF_low.at(i) = CF_low.at(i-1) + frequencies.at(i-1);
    }


    //We also need to know the global cumulative frequency (of all 
    //symbols), which will be the denominator of a formula below.
    //It turns out this value is already stored as CF_low.at(max_symbol+1)
    u64 global_cumulative_frequency = CF_low.at(EOF_SYMBOL+1);
    
    assert(global_cumulative_frequency <= 0xffffffff); //If this fails, frequencies must be scaled down


    u32 lower_bound = 0;
    u32 upper_bound = ~0;

    u32 encoded_bits = 0;

    unsigned int bits_read = 0;

    for(int i = 0; i < 32; i++){
        encoded_bits = (encoded_bits<<1) | stream.read_bit();
        bits_read++;
    }

    
    while(frame_data.size() != frame_len){
        //For safety, we will use u64 for all of our intermediate calculations.
        u64 current_range = (u64)upper_bound - (u64)lower_bound + 1;

        //Figure out which symbol comes next (we can definitely do better than linear 
        //search for this...)

        //First scale the encoded bitstring (which lies between lower_bound and upper_bound)
        //to the range [0, global_cumulative_frequency)
        //With pure real arithmetic, this is equivalent to the equation
        //  scaled = (encoded-low)*(global_cumulative_frequency/current_range),
        //however, we have to salt it with +1 and -1 terms (and rearrange it) to accommodate
        //fixed-point arithmetic.
        u64 scaled_symbol = (((u64)encoded_bits - lower_bound + 1)*global_cumulative_frequency - 1)/current_range;

        u32 symbol = 0;
        while(CF_low.at(symbol+1) <= scaled_symbol)
            symbol++;
            
        //Output the symbol
        frame_data.push_back(symbol - 127);  // - 2047

        //Now that we know what symbol comes next, we repeat the same process as the compressor
        //to prepare for the next iteration.

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

        //Even though we don't have to output bits, we do have to 
        //adjust the lower and upper bounds just like the compressor does.
        while(1){
            //Check if most significant bits (bit index 31) match.
            if ((upper_bound>>31) == (lower_bound>>31)){ 

                //Shift out the MSB of the lower bound, the upper bound and the encoded string
                //(Note that if lower and upper bounds have the same MSB, so does the encoded
                // bitstring)


                //Shift out the MSB of upper_bound (and shift in a 1 from the right)
                upper_bound <<= 1;
                upper_bound |= 1;
                
                //Shift out the MSB of lower_bound (and allow a 0 to be shifted in from the right)
                lower_bound <<= 1;
                
                //Shift out the MSB of encoded_bits (and bring in a new encoded bit from the
                //output file on the right)
                encoded_bits <<= 1;
                if(bits_read == encoded_len) {
                    encoded_bits |= 1;
                }
                else {
                    encoded_bits |= stream.read_bit();
                    bits_read++;
                }
                

            }else if ( ((lower_bound>>30)&0x1) == 1 && ((upper_bound>>30)&0x1) == 0){
                //If the MSBs didn't match, then the MSB of upper_bound must be 1 and
                //the MSB of lower_bound must be 0.
                //If we discover that lower_bound = 01... and upper_bound = 10... 
                //(which is what the if-statement above tests), then we have
                //to account for underflow.

                //If upper_bound = 10(xyz...), set upper_bound = 1(xyz...)
                //(that is, splice out the second-most-significant bit)
                upper_bound <<= 1;
                upper_bound |= (1U<<31);
                upper_bound |= 1;

                //If lower_bound = 01(abc...), set lower_bound = 0(abd...)
                lower_bound <<= 1;
                lower_bound &= (1U<<31) - 1; //i.e. 0x7fffffff

                //Since upper = 10... and lower = 01..., we know that 
                //either encoded_bits = 10... or encoded_bits = 01...
                //(since encoded_bits must be between lower and upper)
                //We want to splice out the second-most-significant bit
                //of encoded_bits (and bring in a new bit on the right)

                //Long way:
                {
                    u32 msb = encoded_bits>>31;
                    u32 rest = encoded_bits&0x3fffffff; //Bits 0 - 30
                    if(bits_read == encoded_len) {
                        encoded_bits = (msb<<31)|(rest<<1)|1;
                    }
                    else {
                        encoded_bits = (msb<<31)|(rest<<1)|stream.read_bit();
                        bits_read++;
                    }
                    
                }
                //Short way (tricky):
                //encoded_bits <<= 1; //Shift everything left (eliminating MSB)
                //encoded_bits = encoded_bits^(1<<31); //Flip bit 31 (which was previously bit 30, which we know was the opposite of the old bit 31)
                //encoded_bits |= stream.read_bit();
            }else{
                break;
            }
        }
    }
}


// lossless compression stuff ends


int main(int argc, char** argv){

    //Note: This program must not take any command line arguments. (Anything
    //      it needs to know about the data must be encoded into the bitstream)
    
    InputBitStream input_stream {std::cin};

    unsigned int quality_code = input_stream.read_u16();
    u32 height {input_stream.read_u32()};
    u32 width {input_stream.read_u32()};

    // quality settings
    if(quality_code == 1) {
        for(unsigned int i = 0; i < 8; ++i){
            for(unsigned int j = 0; j < 8; ++j) {
                Lumi.at(i).at(j) = Lumi.at(i).at(j) * 2;
                Chromi.at(i).at(j) = Chromi.at(i).at(j) * 2;
            }
        }

    }else if(quality_code == 3) {
        for(unsigned int i = 0; i < 8; ++i){
            for(unsigned int j = 0; j < 8; ++j) {
                Lumi.at(i).at(j) = Lumi.at(i).at(j) / 2;
                Chromi.at(i).at(j) = Chromi.at(i).at(j) / 2;
            }
        }
    }
    else if(quality_code == 2) {
        for(unsigned int i = 0; i < 8; ++i){
            for(unsigned int j = 0; j < 8; ++j) {
                Lumi.at(i).at(j) = Lumi.at(i).at(j);
                Chromi.at(i).at(j) = Chromi.at(i).at(j);
            }
        }
    }
    make_C();

    auto Y_prev = create_2d_vector<int>(height,width);
    auto Cb_prev = create_2d_vector<int>(height/2,width/2);
    auto Cr_prev = create_2d_vector<int>(height/2,width/2);

    YUVStreamWriter writer {std::cout, width, height};

    while (input_stream.read_byte()){
        YUVFrame420& frame = writer.frame();

        auto Y_out = create_2d_vector<int>(height,width);
        auto Cr_out = create_2d_vector<int>(height/2, width/2);
        auto Cb_out = create_2d_vector<int>(height/2, width/2);


        unsigned int frame_type = input_stream.read_bit();
        unsigned int motion_vector_count = input_stream.read_u32();
        
        std::vector<std::vector<unsigned int>> motion_vectors {};
        unsigned int motion_vec_ind = 0;
        for(unsigned int i = 0; i < motion_vector_count; ++i) {
            std::vector<unsigned int> temp {};
            temp.push_back(input_stream.read_u32());
            temp.push_back(input_stream.read_u32());
            motion_vectors.push_back(temp);
        }
        
        unsigned int frame_len = input_stream.read_u32();
        unsigned int encoded_len = input_stream.read_u32();

        std::vector<int> frame_data;
        arith_decompress(input_stream, frame_data, frame_len, encoded_len);
        

        int frame_curr_ind = 0;
        frame_data = rle_decode(frame_data);

        // Pre_block input:
        //  I-Frame: 256 Y elements -> 64 Cr elements -> 64 Cb elements
        //  P-Frame: 2 values for macroblock motion vector -> 256 Y elements -> 64 Cr elements -> 64 Cb elements
        for(unsigned int row = 0; row < height; row+=16) {
            for (unsigned int col = 0; col < width; col+=16) {
                
                std::vector<unsigned int> local_motion_vector(2, 0);
                if(frame_type == 1) {
                    local_motion_vector = motion_vectors.at(motion_vec_ind);
                    motion_vec_ind++;  
                }
                
                // read a yuv 4:2:0 macroblocks
                //std::vector<std::vector<int>> Y_blocks(4, std::vector<int>());
                auto Y_block = create_2d_vector<int>(16,16);
                auto Cb_block = create_2d_vector<int>(8,8);
                auto Cr_block = create_2d_vector<int>(8,8);

                // rebuild 16x16 Y block
                for(unsigned int r = 0; r < 16; r+=8) {
                    for(unsigned int c = 0; c < 16; c+=8) {
                        for(unsigned int i = 0; i < 64; ++i) {
                            Y_block.at(r + Zigzag.at(i).first).at(c + Zigzag.at(i).second) = frame_data.at(frame_curr_ind);
                            frame_curr_ind++;
                        }
                    }
                }
                // rebuild 8x8 Cb block
                for(unsigned int i = 0; i < 64; ++i) {
                    Cb_block.at(Zigzag.at(i).first).at(Zigzag.at(i).second) = frame_data.at(frame_curr_ind);
                    frame_curr_ind++;
                }
                // rebuild 8x8 Cr block
                for(unsigned int i = 0; i < 64; ++i) {
                    Cr_block.at(Zigzag.at(i).first).at(Zigzag.at(i).second) = frame_data.at(frame_curr_ind);
                    frame_curr_ind++;
                }

                rev_block_DCT(Y_block, false);
                rev_block_DCT(Cb_block, true);
                rev_block_DCT(Cr_block, true);

                if(frame_type == 1) {   // P-Frame
                    // format block back into the proper location in frame and add prev frame values at motion vector
                    place_P_frame_block(Y_block, Y_prev, Y_out, frame, local_motion_vector.at(0), local_motion_vector.at(1), height, width, row, col, 0);
                    place_P_frame_block(Cb_block, Cb_prev, Cb_out, frame, local_motion_vector.at(0)/2, local_motion_vector.at(1)/2, height/2, width/2, row/2, col/2, 1);
                    place_P_frame_block(Cr_block, Cr_prev, Cr_out, frame, local_motion_vector.at(0)/2, local_motion_vector.at(1)/2, height/2, width/2, row/2, col/2, 2);
                    
                }
                else {  // I-Frame
                    // format block back into the proper location in frame
                    place_I_frame_block(Y_block, Y_out, frame, height, width, row, col, 0);
                    place_I_frame_block(Cb_block, Cb_out, frame, height/2, width/2, row/2, col/2, 1);
                    place_I_frame_block(Cr_block, Cr_out, frame, height/2, width/2, row/2, col/2, 2);
                }
            }
        }

        Y_prev = Y_out;
        Cb_prev = Cb_out;
        Cr_prev = Cr_out;

        // output frame
        writer.write_frame();
 
    }


    return 0;
}