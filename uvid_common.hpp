/* uvid_common.hpp
   
   Jonathan Cote, V00962634
   
*/
#ifndef uvid_common_hpp
#define uvid_common_hpp
#include <vector>
#include <string>
#include <math.h>


std::vector<std::pair<unsigned int, unsigned int>> Zigzag = 
    {
        {0,0}, 
        {0,1}, {1,0}, 
        {2,0}, {1,1}, {0,2}, 
        {0,3}, {1,2}, {2,1}, {3,0}, 
        {4,0}, {3,1}, {2,2}, {1,3}, {0,4}, 
        {0,5}, {1,4}, {2,3}, {3,2}, {4,1}, {5,0},
        {6,0}, {5,1}, {4,2}, {3,3}, {2,4}, {1,5}, {0,6},
        {0,7}, {1,6}, {2,5}, {3,4}, {4,3}, {5,2}, {6,1}, {7,0},
        {7,1}, {6,2}, {5,3}, {4,4}, {3,5}, {2,6}, {1,7},
        {2,7}, {3,6}, {4,5}, {5,4}, {6,3}, {7,2},
        {7,3}, {6,4}, {5,5}, {4,6}, {3,7},
        {4,7}, {5,6}, {6,5}, {7,4},
        {7,5}, {6,6}, {5,7},
        {6,7}, {7,6},
        {7,7}
    };



// provide from A3 uvg_common by B. Bird
template<typename T>
std::vector<std::vector<T> > create_2d_vector(unsigned int outer, unsigned int inner){
    std::vector<std::vector<T> > V {outer, std::vector<T>(inner,T() )};
    return V;
}

// provide from A3 uvg_common by B. Bird
//The floating point calculations we use while converting between 
//RGB and YCbCr can occasionally yield values slightly out of range
//for an unsigned char (e.g. -1 or 255.9).
//Furthermore, we want to ensure that any conversion uses rounding
//and not truncation (to improve accuracy).
inline unsigned char round_and_clamp_to_char(double v){
    //Round to int 
    int i = (int)(v+0.5);
    //Clamp to the range [0,255]
    if (i < 0)
        return 0;
    else if (i > 255)
        return 255;
    return i;
}


void block_zigzag(std::vector<std::vector<int>>& in, std::vector<int>& zigzag_block) {
    // Convert to ZigZag storage pattern
    for(unsigned int i = 0; i < 64; ++ i) {
        zigzag_block.push_back(in.at(Zigzag.at(i).first).at(Zigzag.at(i).second));
    }
}

void block_de_zigzag(std::vector<int>& block, std::vector<std::vector<int>>& out) {
    // de-zig-zag
    for(unsigned int i = 0; i < 64; ++i) {
        out.at(Zigzag.at(i).first).at(Zigzag.at(i).second) = block.at(i);
    }
}


std::vector<int> rle_encode(std::vector<int>& stream) {
    // simple RLE
    std::vector<int> out {};
    int count = 0;
    int curr;
    
    for(unsigned int i = 0; i < stream.size(); ++i) {
        curr = stream.at(i);
        if(curr == 0) {
            if(count > 125) {
                out.push_back(0);
                out.push_back(count);
                count = 0;
            }
            count++;
        }
        else if(count > 0) {
            out.push_back(0);
            out.push_back(count);
            count = 0;
            out.push_back(curr);
        }
        else {
            out.push_back(curr);
        }
    }

    if(count > 0) {
        out.push_back(0);
        out.push_back(count);
    }

    return out;

}

std::vector<int> rle_decode(std::vector<int>& stream) {
    // simple RLE
    std::vector<int> out {};
    //int count = 0;
    bool zero_run = false;

    for(unsigned int i = 0; i < stream.size(); ++i) {
        if(zero_run) {
            for(int j = 0; j < stream.at(i); ++j) {
                out.push_back(0);
            }
            zero_run = false;
        }
        else if(stream.at(i) == 0) {
            zero_run = true;
        }
        else {
            out.push_back(stream.at(i));
        }
    }
   
    return out;

}


std::vector<std::vector<int>> get_block(std::vector<std::vector<int>>& frame, unsigned int height,
                                        unsigned int width, unsigned int dx, unsigned int dy, unsigned int block_size) 
{
    auto out_block = create_2d_vector<int>(block_size, block_size);
    for(unsigned int y = 0; y < block_size; ++y) {
        for(unsigned int x = 0; x < block_size; ++x) {
            unsigned int frame_row = y + dy;
            unsigned int frame_col = x + dx;

            if(frame_row >= height) {
                frame_row = height - 1;  // Copy last row
            }
            if(frame_col >= width) {
                frame_col = width - 1;  // Copy last col
            }

            out_block.at(y).at(x) = frame.at(frame_row).at(frame_col);
        }
    }
    return out_block;
}


void revs_Block_Quantized_DCT(std::vector<std::vector<int>>& block_in, std::vector<std::vector<int>>& block_out,
                              std::vector<std::vector<int>>& quantizer, std::vector<std::vector<double>>& C,
                              std::vector<std::vector<double>>& C_T, unsigned int x_start, unsigned int y_start)
{
    auto sub_De_Quant = create_2d_vector<double>(8,8);
    auto sub_temp = create_2d_vector<double>(8,8);

    // reverse Quantization
    for(unsigned int i = 0; i < 8; ++i){
        for(unsigned int j = 0; j < 8; ++j) {
            sub_De_Quant.at(i).at(j) = static_cast<double>(block_in.at(y_start + i).at(x_start + j) * quantizer.at(i).at(j));
        } 
    }

    // reverse DCT
    for(unsigned int i = 0; i < 8; ++i) {
        for(unsigned int j = 0; j < 8; ++j) {
            double sum = 0;
            for(unsigned int k = 0; k < 8; ++k) {
                sum += C_T.at(i).at(k) * sub_De_Quant.at(k).at(j);
            }
            sub_temp.at(i).at(j) = sum;
        }
    } 

    for(unsigned int i = 0; i < 8; ++i) {
        for(unsigned int j = 0; j < 8; ++j) {
            double sum = 0;
            for(unsigned int k = 0; k < 8; ++k) {
                sum += sub_temp.at(i).at(k) * C.at(k).at(j);
            }
            block_out.at(y_start + i).at(x_start + j) = static_cast<int>(sum + 0.5);
        }
    } 
}


void calc_Block_Quantized_DCT(std::vector<std::vector<int>>& block_in, std::vector<std::vector<int>>& block_out,
                              std::vector<std::vector<int>>& quantizer, std::vector<std::vector<double>>& C,
                              std::vector<std::vector<double>>& C_T, unsigned int x_start, unsigned int y_start)
{  
    
    auto C_sub = create_2d_vector<double>(8,8);
    auto sub_DCT = create_2d_vector<double>(8,8);

    // DCT matrix multiplication
    for(unsigned int i = 0; i < 8; ++i) {
        for(unsigned int j = 0; j < 8; ++j) {
            double sum = 0;
            for(unsigned int k = 0; k < 8; ++k) {
                sum += C.at(i).at(k) * block_in.at(y_start + k).at(x_start + j);
            }
            C_sub.at(i).at(j) = sum;
        }
    }
    for(unsigned int i = 0; i < 8; ++i) {
        for(unsigned int j = 0; j < 8; ++j) {
            double sum = 0;
            for(unsigned int k = 0; k < 8; ++k) {
                sum += C_sub.at(i).at(k) * C_T.at(k).at(j);
            }
            sub_DCT.at(i).at(j) = sum;
        }
    }

    // Quantization
    for(unsigned int i = 0; i < 8; ++i){
        for(unsigned int j = 0; j < 8; ++j) {
            int q_val = std::max(std::min(static_cast<int>((sub_DCT.at(i).at(j) / quantizer.at(i).at(j)) + 0.5), 127), -127);
            block_out.at(y_start + i).at(x_start + j) = q_val;
        }     
    }
}





#endif