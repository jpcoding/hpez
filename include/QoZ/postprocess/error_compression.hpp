#ifndef SZ_ERROR_COMPRESSION_HPP
#define SZ_ERROR_COMPRESSION_HPP
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include<vector> 
#include<cmath>
#include<limits>
#include <memory>
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/utils/ByteUtil.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/utils/FileUtil.hpp"
#include "QoZ/utils/Interpolators.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/utils/Timer.hpp"
#include "QoZ/quantizer/Quantizer.hpp"
#include "QoZ/encoder/Encoder.hpp"
#include "QoZ/lossless/Lossless.hpp"
#include "QoZ/encoder/Encoder.hpp"
#include "QoZ/encoder/HuffmanEncoder.hpp"
#include "QoZ/lossless/Lossless.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/predictor/Predictor.hpp"
#include "QoZ/utils/Interpolators.hpp"
// #include "QoZ/compressor/SZInterpolationCompressor.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
#include "QoZ/postprocess/error_interp_compressor.hpp"



namespace QoZ {


template <class T, uint N>



class ErrorCompressor {


public: 
  ErrorCompressor(size_t orig_num_elements, const size_t* global_dimensions,
    const size_t* dimension_offsets, const int interp_id)
    {
        for(int i = 0; i < N; i++)
        {
            orig_global_dimensions[i] = global_dimensions[i];
            orig_dimension_offsets[i] = dimension_offsets[i];
        }
        interpolation_order_id = interp_id;
        

    }



void interpolate_error_slice(T* error_smaple_begin, T* interpolated_slice, int sample_stride,
                int interpolated_dim1, int interpolated_dim2, 
                            int plane_dir1, int plane_dir2)

{
    // dim1 is the fast changing dimension 
    // using bilinear interpolation 
    int sample_dim1 = error_global_dimensions[plane_dir1];
    int sample_dim2 = error_global_dimensions[plane_dir2];
    for(int i = 0; i < sample_dim1-1; i+=1)
    {
        for(int j = 0; j < sample_dim2-1; j+=1)
        {
            // get the four corner points 
            T p00 = error_smaple_begin[i * error_dimension_offsets[plane_dir1] + j *error_dimension_offsets[plane_dir2] ]; // x1, y1
            T p01 = error_smaple_begin[i * error_dimension_offsets[plane_dir1] + (j+1) *error_dimension_offsets[plane_dir2] ]; // x1, y2
            T p10 = error_smaple_begin[(i+1) * error_dimension_offsets[plane_dir1] + (j) *error_dimension_offsets[plane_dir2] ]; // x2, y1
            T p11 = error_smaple_begin[(i+1) * error_dimension_offsets[plane_dir1] + (j+1) *error_dimension_offsets[plane_dir2] ]; // x2, y2
            for(int ii = 0; ii < sample_stride; ii++)
            {
                for(int jj = 0; jj < sample_stride; jj++)
                {
                    int plane_index = (i * sample_stride + ii) *interpolated_dim2  + j * sample_stride + jj;
                    interpolated_slice[plane_index] = 1.0/(sample_stride*sample_stride) * (
                        p00 * (sample_stride - ii) * (sample_stride - jj) +
                        p01 * (sample_stride - ii) * jj + 
                        p10 * ii * (sample_stride - jj) + 
                        p11 * ii * jj);
                }
            }
  
        }
    }
}

void compensate_with_interpolation(T* dec_data, T *dec_error_downsampled, int slice_start_idx,
                int interp_direcion, int interp_dir_stride, 
                int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_sample_stride, int interp_dir_sample_stride=1)
{
    // 1. get downsampled 



}

void compensate_with_interpolated_slice(T* dec_data_begin, T* interpolated_slice,
                int interp_direcion, int interp_dir_stride, 
                int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_dir1_dim, int plane_dir2_dim)
{
    // get the start index of the slice 
    // size_t current_global_index = slice_idx * interp_dir_stride * interp_dir_sample_stride * orig_global_dimensions[interp_direcion];
    T* dec_data_ptr = dec_data_begin;
    T* interp_slice_ptr = interpolated_slice;
    // std::cout<< "plane_dir1_dim: " << plane_dir1_dim << std::endl;
    // std::cout << "plane_dir2_dim: " << plane_dir2_dim << std::endl;

    // std::cout << "plane global dim 1: " << orig_global_dimensions[plane_dir1] << std::endl;
    // std::cout << "plane global dim 2: " << orig_global_dimensions[plane_dir2] << std::endl;
    for(int i = 0; i < plane_dir1_dim; i++)
    {
        for(int j = 0; j < plane_dir2_dim; j++)
        {
           dec_data_ptr = dec_data_begin + i * plane_dir1_stride * orig_dimension_offsets[plane_dir1] 
                            + j * plane_dir2_stride * orig_dimension_offsets[plane_dir2];
           interp_slice_ptr = interpolated_slice + i * plane_dir2_dim + j;
              *dec_data_ptr += *interp_slice_ptr;
        }
    }
}






std::vector<T> downsample_error_3d(const T *orig_data, T*dec_data, 
                int interp_direcion, int interp_dir_stride, 
                int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_sample_stride, int interp_dir_sample_stride=1) {
    // define: error = orig - dec 
    // assert N = 3 
    static_assert(N == 3, "N must be 3");
    // prepare the downsampled error dimensions and offsets 
    error_global_dimensions[plane_dir1] = orig_global_dimensions[plane_dir1]/plane_dir1_stride/plane_sample_stride;
    error_global_dimensions[plane_dir2] = orig_global_dimensions[plane_dir2]/plane_dir2_stride/plane_sample_stride;
    error_global_dimensions[interp_direcion] = orig_global_dimensions[interp_direcion]/interp_dir_stride/interp_dir_sample_stride;

    // offsets 
    error_dimension_offsets[2] = 1;
    error_dimension_offsets[1] = error_global_dimensions[2];
    error_dimension_offsets[0] = error_global_dimensions[1] * error_global_dimensions[2];
    
    // std::cout << "error_global_dimensions: " << error_global_dimensions[0] << " " << error_global_dimensions[1] << " " << error_global_dimensions[2] << std::endl;

    // init the error sample storage
    size_t error_sample_size = 1;
    for (int i = 0; i < N; i++) {
        error_sample_size *= error_global_dimensions[i];
    }
    std::vector<T> error_sample(error_sample_size, 0); 

    // downsample error along the interpolation direction 
    size_t current_global_index, current_error_index;
    for(int i = 0; i < error_global_dimensions[interp_direcion]; i+= 1) {
        for(int j = 0; j < error_global_dimensions[plane_dir1]; j++) {
            for(int k = 0; k < error_global_dimensions[plane_dir2]; k++) {
                current_global_index = (i*interp_dir_stride +1)*interp_dir_sample_stride*orig_dimension_offsets[interp_direcion] +
                j*plane_dir1_stride*plane_sample_stride*orig_dimension_offsets[plane_dir1] + k*plane_dir2_stride*plane_sample_stride*orig_dimension_offsets[plane_dir2];
                current_error_index = i*error_dimension_offsets[interp_direcion] + j*error_dimension_offsets[plane_dir1] + k*error_dimension_offsets[plane_dir2];
                error_sample[current_error_index] = orig_data[current_global_index] - dec_data[current_global_index];
                // printf("global x %d, y %d, z %d\n", (i*interp_dir_stride+1), j*plane_dir1_stride, k*plane_dir2_stride);
                // printf("error x %d, y %d, z %d\n", i, j, k);
                // error_sample[current_error_index] = orig_data[current_global_index];
            }
        }
    }
    // std::cout <<"max error " << *std::max_element(error_sample.begin(), error_sample.end()) << std::endl;
    // std::cout <<"min error " << *std::min_element(error_sample.begin(), error_sample.end()) << std::endl;

    return error_sample; 
}




/*
interpolation direciton is perperdicular to the plane slices
sample stride if for the slice
// */
// uchar *compress_error_3d(const T *orig_data, T*dec_data,  size_t &compressed_error_size , int interp_direcion, int interp_dir_stride,int plane_dir1, int plane_dir2,
//                 int plane_dir1_stride, int plane_dir2_stride,
//                 int plane_sample_stride, int interp_dir_sample_stride=1) {
    // 
    // auto error_sample = downsample_error_3d(orig_data, dec_data, interp_direcion, interp_dir_stride, plane_dir1, plane_dir2, plane_dir1_stride, plane_dir2_stride, plane_sample_stride, interp_dir_sample_stride);
    // double rel_eb = 0.005;
    // double abs_eb = (*std::max_element(error_sample.begin(), error_sample.end()) - *std::min_element(error_sample.begin(), error_sample.end()))*rel_eb;
    // // build a quantizer
    // // auto error_quantizer = QoZ::LinearQuantizer<T>(abs_eb);
    // // auto error_encoder = QoZ::HuffmanEncoder<int>();
    // // auto error_lossless = QoZ::Lossless_zstd();

    // auto error_config = QoZ::Config(error_global_dimensions[0], error_global_dimensions[1], error_global_dimensions[2]);
    // error_config.error_smoothing = false;
    // error_config.compress_error = false;
    // error_config.absErrorBound = abs_eb; 

    // auto interpolation_compressor = SZErrorInterpolationCompressor<T, 3, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
    //         QoZ::LinearQuantizer<T>(error_config.absErrorBound, error_config.quantbinCnt / 2),
    //         QoZ::HuffmanEncoder<int>(),
    //         QoZ::Lossless_zstd());

    // // compress the error sample 
    // // build a compressor 
    // // auto interpolation_compressor = QoZ::SZInterpolationCompressor<T,3,QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd> (error_quantizer, error_encoder,error_lossless);
    // // make a config file 


    // compressed_error_size = 0;
    // uchar* compressed_error = interpolation_compressor.compress(error_config, error_sample.data(), compressed_error_size, false);

//     uchar* compressed_error = nullptr;

//     return 0;
// }

uchar* compress_error_3d(T* downsampled_error, size_t &compressed_error_size, double abs_eb )
{
    // auto error_config = QoZ::Conclearfig(error_global_dimensions[0], error_global_dimensions[1], error_global_dimensions[2]);
    // error_config.error_smoothing = false;
    // error_config.compress_error = false;
    // error_config.absErrorBound = abs_eb; 
    // error_config.error_smoothing = false;
    // error_config.compress_error = false;
    // error_config.quantization_prediction_on = 0;
    // error_config.quantization_prediction_start_level= 0;
    // error_config.blockSize = 9999;
    // // std::cout << "error_config.absErrorBound = " << error_config.absErrorBound << std::endl; 

    // // compress the error sample 
    // // build a compressor 
    // std::cout << "error_config.num = " << error_config.num << std::endl;
    
    // auto interpolation_compressor = SZErrorInterpolationCompressor<T, 3, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
    //         QoZ::LinearQuantizer<T>(error_config.absErrorBound, error_config.quantbinCnt),
    //         QoZ::HuffmanEncoder<int>(),
    //         QoZ::Lossless_zstd());
    //             // make a config file 
    

    // compressed_error_size = 0;
    // uchar* compressed_error = interpolation_compressor.compress(error_config, downsampled_error, compressed_error_size);
    //     // uchar *mmmm = new uchar[60];   
    // return compressed_error;
    return nullptr; 
}

template <uint NN = N>
typename std::enable_if<NN == 1, uchar*>::type
error_compensation_3d(const T *orig_data, T*dec_data,  size_t &compressed_error_size ,
                int interp_direcion, int interp_dir_stride,int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_sample_stride, int interp_dir_sample_stride=1)
{
    return nullptr;

}

template <uint NN = N>
typename std::enable_if<NN == 2, uchar*>::type
error_compensation_3d(const T *orig_data, T*dec_data,  size_t &compressed_error_size ,
                int interp_direcion, int interp_dir_stride,int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_sample_stride, int interp_dir_sample_stride=1)
{

return nullptr;
}

template <uint NN = N>
typename std::enable_if<NN == 4, uchar*>::type
error_compensation_3d(const T *orig_data, T*dec_data,  size_t &compressed_error_size ,
                int interp_direcion, int interp_dir_stride,int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_sample_stride, int interp_dir_sample_stride=1)
{
return nullptr;
}




template <uint NN = N>
typename std::enable_if<NN == 3, uchar*>::type
error_compensation_3d(const T *orig_data, T*dec_data,  size_t &compressed_error_size ,
                int interp_direcion, int interp_dir_stride,int plane_dir1, int plane_dir2,
                int plane_dir1_stride, int plane_dir2_stride,
                int plane_sample_stride, int interp_dir_sample_stride=1)
{
        //1. get sampled  error
        auto error_sample = downsample_error_3d(orig_data, dec_data, interp_direcion, interp_dir_stride, 
        plane_dir1, plane_dir2, plane_dir1_stride, plane_dir2_stride, plane_sample_stride, interp_dir_sample_stride);
        //2. compress the error
        compressed_error_size = 0; 
        double rel_eb = 0.01;
        double abs_eb = (*std::max_element(error_sample.begin(), 
                        error_sample.end()) - *std::min_element(error_sample.begin(), error_sample.end()))*rel_eb;
        // std::cout<< "abs_eb = " << abs_eb << std::endl;
        uchar* compressed_error = compress_error_3d(error_sample.data(),
                                                    compressed_error_size,abs_eb );

        // std::cout << "compressed_error_size = "<< compressed_error_size  << std::endl;


        //3. interpolate the error slice by slice and compensate the dec_data
        // 3.1 create a slice of the original size 
        int plane_dir1_dim = orig_global_dimensions[plane_dir1]/plane_dir1_stride;
        int plane_dir2_dim = orig_global_dimensions[plane_dir2]/plane_dir2_stride;
        std::vector<T> interpolated_slice(plane_dir1_dim*plane_dir2_dim, 0);
        int i ; 
        for( i = 0; i<error_global_dimensions[interp_direcion]; i++)
        {
            std::fill(interpolated_slice.begin(), interpolated_slice.end(), 0);
            T* error_sample_slice_begin = error_sample.data() + i * error_dimension_offsets[interp_direcion];
            // // interpolate the slice 
            interpolate_error_slice( error_sample_slice_begin, interpolated_slice.data(),plane_sample_stride,
                                            plane_dir1_dim, plane_dir2_dim, plane_dir1, plane_dir2);
            // compensate the dec_data
            T* dec_data_begin = dec_data + (i * interp_dir_stride+1) * interp_dir_sample_stride * orig_dimension_offsets[interp_direcion]; 
            compensate_with_interpolated_slice(dec_data_begin, interpolated_slice.data(),
                                                interp_direcion, interp_dir_stride, plane_dir1, plane_dir2, 
                                                plane_dir1_stride, plane_dir2_stride, plane_dir1_dim, plane_dir2_dim);
        }
        // std::cout << "i = " << i << std::endl;

        // std::cout << compressed_error << std::endl; 
        // std::cout << "completed error compensation" << std::endl;
        return compressed_error;
}



private:


void init()
{
    for(int i = 0; i < N; i++)
    {
        error_global_dimensions[i] = orig_global_dimensions[i];
        error_dimension_offsets[i] = orig_dimension_offsets[i];
    }

    // auto interpolation_compressor = QoZ::SZInterpolationCompressor(quantizer, encoder, lossless);


}
// int sample_stride;
// size_t error_sample_size;
// std::vector<T> error_sample;


size_t orig_num_elements;
std::array<size_t, N> orig_global_dimensions; // selected from 
std::array<size_t, N> orig_dimension_offsets;
std::vector<std::array<int, N>> orig_dimension_sequences;
int interpolation_order_id; 
std::array<size_t, N> orig_selected_dimensions; 

std::array<size_t, N> error_global_dimensions;
std::array<size_t, N> error_dimension_offsets;
std::vector<std::array<int, N>> error_dimension_sequences;

std::array<size_t,N> interpolated_error_global_dimensions;
std::array<size_t,N> interpolated_error_dimension_offsets;
std::vector<std::array<int, N>> interpolated_error_dimension_sequences;

//compression utils 
// Quantizer quantizer;
// Encoder encoder;
// Lossless lossless;


};


};


#endif // SZ_ERROR_COMPRESSION_HPP