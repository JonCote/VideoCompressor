Jonathan Cote V00962634


Overview:

    Video compression frame by frame using 8x8 blocking DCT quantization. Compressor loops through each frame of the video compressing the frame one marcroblock (four 8x8 Y channel blocks, one 8x8 Cb block, and one 8x8 Cr block) at a time. After aquiring a set of blocks for the YUV macroblock the compressor checks if current frame is a I-frame or P-frame. If current frame is a I-frame the YUV macroblocks will go through DCT, quantization, and block reordering using zigzaging before being stored for writing to the output file once frame compression is complete. If the current frame is a P-frame the compressor determines a motion vector for the YUV macroblock using a binary/logarithmic search. Once a suitable motion vector is found the compressor takes the YUV macroblocks and corresponding previous frames blocks found using the motion vector to determines the current block difference from corresponding previous frame block. The compressor then takes the difference blocks and runs them through the same DCT, quantization, and zigzaging process as a I-Frame. Once ever YUV macroblock is processed and stored for the frame a simple RLE processes is run on the frames data to compress runs of 0's. Finally, the compressor runs the frames lossy compressed data through a adaptive arithmetic compression before pushing the compressed frame data to the output file. The above process is repeated for each frame of the video where the first and ever 10th frame will be compressed as a I-frame and all other frames as P-frames.

    Decompression follows the same pipeline in reverse with out motion vector search and determines frame type by reading a flag value sent before each frame in the compressed data stream.

    compression tests: (ran on remote to class server)
        sintel_trailer_2k_720p24.y4m (1280x720 24fps 1253frames 1,691,558KB) compressed in 6mins to 55,727KB (compression ratio of 30x)
        sintel_trailer_2k_480p24.y4m (854x480 24fps 1253frames 752,395KB) compressed in 3mins to 26,442KB (compression ratio of 28x)
        deadline_cif.y4m (352x288 30fps 1374frames 204,048KB) compressed in 1min to 9,714KB (compression ratio of 21x)


Features Implemented:

    - basic requirements
    - motion vector calculation per macroblock of P-frame
    - compressed bitstream using a dynamic entropy coding technique (adaptive arithmetic encoding)
    - implemented an advanced search technique for motion vectors (binary/logarithmic search more detail below)
    - significantly higher compression ratio than required (did not talk to instructor so if not valid it is ok)



Architecture:

    DCT quantization video compression pipeline using I-frames and P-frames. compression pipeline is ran frame by frame processing one YUV macroblock consisting of four 8x8 Y channel blocks, one 8x8 Cb channel block, and one 8x8 Cr channel block at a time. Each frame of the video is marked as a I-frame or P-frame, compression pipeline is dependent on frame type and is expressed in more detail below.

    Components/structures:
        - YUV macroblock collection: grabs the next set of 8x8 blocks to be processed for the frame storing each block into a vector
        
        - DCT + Quantization: 
            - 2d-DCT on 8x8 block followed by Quantization using Quantization table specified by quality setting
            - reverse DCT + Quantization stored in previous frame data, done during compression to limit errors in P-frames during decompression

        - Block zigzag reordering:
            - zigzag reordering of block data to improve data locality for RLE compression
        
        - Binary/logarithmic motion vector search: 
            - search previous frame for best matching block to current block.
            - before binary search is conducted a min_AAD value is calculated using the previous block at same location and returns this location if AAD is <= 5.0. 
            - binary search is conducted by subdividing the previous frame into 4 smaller sub-frames and calculating the average absolute difference (AAD) of the center 8x8 block in the sub-frame with the current 8x8 block. The sub-frame with the best AAD is kept then subdivided. subdivision is continued until subdivision can not produce min 8x8 sub-frame sizes.
            - binary search will return the motion vector with the best AAD

        - RLE encoding: 
            - simple RLE encoding on 0 value runs, with runlength cap of 125
            - encoded runs as 0 followed by count
            - decoding if 0 is seen next value tells the number of zeros it expands to
        
        - Arithmetic coding: 
            - adaptive arithmetic encoding
            - updates frequency table as it encodes/decodes each character
            - used arithmetic encoder provided by B. Bird with adaption to support frequency value updating and frame by frame block encoding/decoding
        
    
    I-Frame pipeline:
        Per YUV macroblock: 
            collect 8x8 blocks -> DCT + Quantization -> reverse DCT + Quantization (for previous frame) -> block zigzag reordering -> store in frame_out
    
    P-Frame pipeline:
        Per YUV macroblock:
            collect 8x8 blocks -> per block motion vector binary search -> block difference calculation -> DCT + Quantization -> reverse DCT + Quantization (for previous frame) -> block zigzag reordering -> store in frame_out

    Frame data output:
        Once all YUV macroblock for the frame are processed:
            update previous frame data -> RLE encode -> arithmetic coding -> output encoded frame data


    Decompressor: 
        - follows same processes as compressor but in reverse (excluding motion vector search).


Bitstream:

    quality code: u16 value
    height: u32 value
    width: u32 value

    Per frame:
        frame available flag: 
            if another frame: byte representation of integer 1
            else: byte representation of integer 0 (end of file flag)
        
        frame type flag:
            if I-Frame: bit 0
            if P-Frame: bit 1

        - If P-Frame:
            - u32 value corresponding to the number of motion vectors used for the frame
            - per motion vector:
                - 2 u32 values corresponding to blocks location (dx, dy)

        encoded frame size: u32 value, corresponding to the number of data points to be aquired from compressed data file after undoing arithmetic coding

        data block size: u32 value, corresponding to the bits of data stored for the frame in compressed data file

        data block: bit stream of length data block size that is arithmetic and RLE encoded

        data block (after Arithmetic decoding and RLE decoding of the data block): 
            - arithmetic+RLE encoded block that decodes to the following structure,
                - Data is written one YUV macroblock at a time, for each macroblock in the frame the follow data is written (in exact order):
                    - zigzag ordered macroblock data: (repeats tell all macroblocks for frame are written)
                        - 256 values corresponding to the 4 8x8 Y channel blocks
                        - 64 values corresponding to the Cb channel block
                        - 64 values corresponding to the Cr channel block
                    

    end of file flag: byte representation of integer 0 (same as expressed in frame available flag)
    
    extra bits to aquire even byte


Bibliography:

    - arithmetic coding adapted from provided arithmetic coding file from B. Bird
    - create_2d_vector and round_and_clamp_to_char aquired from A3 uvg_common by B. Bird
   

