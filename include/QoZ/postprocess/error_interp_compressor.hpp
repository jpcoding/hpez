#ifndef _SZ_ERROR_INTERPOLATION_COMPRESSSOR_HPP
#define _SZ_ERROR_INTERPOLATION_COMPRESSSOR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include "QoZ/utils/Statistic.hpp"

// #include "QoZ/compressor/SZInterpolation_postprocess.hpp"
#include "QoZ/postprocess/compensate_2d.hpp"
#include "QoZ/postprocess/post_compensate_zero_quant.hpp"
#include "QoZ/def.hpp"
#include "QoZ/encoder/Encoder.hpp"
#include "QoZ/encoder/HuffmanEncoder.hpp"
#include "QoZ/lossless/Lossless.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/predictor/Predictor.hpp"
#include "QoZ/quantizer/Quantizer.hpp"
// #include "QoZ/utils/Accumulator.hpp"
#include "QoZ/utils/ByteUtil.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/utils/FileUtil.hpp"
#include "QoZ/utils/Interpolators.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/utils/Timer.hpp"


namespace QoZ {
template <class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZErrorInterpolationCompressor {
 public:
  SZErrorInterpolationCompressor(
      Quantizer quantizer, Encoder encoder, Lossless lossless) :
      quantizer(quantizer), encoder(encoder), lossless(lossless)
  {
    static_assert(
        std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
        "must implement the quatizer interface");
    static_assert(
        std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
        "must implement the encoder interface");
    static_assert(
        std::is_base_of<concepts::LosslessInterface, Lossless>::value,
        "must implement the lossless interface");
  }

  T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num)
  {
    T *dec_data = new T[num];
    return decompress(cmpData, cmpSize, dec_data);
  }

  T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData)
  {
    size_t remaining_length = cmpSize;
    uchar *buffer = lossless.decompress(cmpData, remaining_length);
    uchar const *buffer_pos = buffer;

    read(global_dimensions.data(), N, buffer_pos, remaining_length);
    read(blocksize, buffer_pos, remaining_length);
    read(interpolator_id, buffer_pos, remaining_length);
    read(direction_sequence_id, buffer_pos, remaining_length);

    // std::cout<<"num_detection_block = " << num_detection_block <<std::endl;
    Timer timer;
    timer.start();


    read(original_max, buffer_pos);
    read(original_min, buffer_pos);
    // std::cout << "detection_eb_rate = " << detection_eb_rate << std::endl;
    // std::cout << "noise_rate = " << noise_rate << std::endl;

    original_range = original_max - original_min;

    init();

    // timer.stop("read auxilliary data");



    quantizer.load(buffer_pos, remaining_length);
    encoder.load(buffer_pos, remaining_length);
    quant_inds =
        encoder.decode(buffer_pos, num_elements);

    encoder.postprocess_decode();

    lossless.postdecompress_data(buffer);
    double eb_input = quantizer.get_eb();
    double eb_final;
    // double eb_reduction_factor;
    if (interpolator_id == 0) {
      eb_final = eb_input /
                 pow(linear_interp_eb_factor, (interpolation_level - 1) * N);
      eb_reduction_factor = pow(linear_interp_eb_factor, N);
    }
    else {
      eb_final = eb_input /
                 pow(cubic_interp_eb_factor, (interpolation_level - 1) * N);
      eb_reduction_factor = pow(cubic_interp_eb_factor, N);
    }

    std::cout << "start decompression\n";
    // *decData = quantizer.recover(0, quant_inds[quant_index++]);
    recover(0, *decData, 0);



    for (uint level = interpolation_level;
         level > 0 && level <= interpolation_level; level--) {
      // if (level >= 3) {
      //     quantizer.set_eb(eb * eb_ratio);
      // } else {
      //     quantizer.set_eb(eb);
      // }

      // direction_sequence_id = direction_choice(mt);
      // std::cout << "direction_sequence_id = " << direction_sequence_id <<
      // std::endl; for (int i = 0; i < N; i++)
      // {
      //     std::cout << dimension_sequences[direction_sequence_id][i] << " ";
      // }
      // std::cout << std::endl;

      // current_level = level;
      // current_base_eb = eb_final;
      // quantizer.set_eb(eb_final);
      // eb_final *= eb_reduction_factor;



      // if (level <= 2) {
      //   quantizer.set_eb(eb_input);
      // }

      size_t stride = 1U << (level - 1);
      auto inter_block_range =
          std::make_shared<QoZ::multi_dimensional_range<T, N>>(
              decData, std::begin(global_dimensions),
              std::end(global_dimensions), stride * blocksize, 0);
      auto inter_begin = inter_block_range->begin();
      auto inter_end = inter_block_range->end();
      for (auto block = inter_begin; block != inter_end; ++block) {
        auto end_idx = block.get_global_index();
        for (int i = 0; i < N; i++) {
          end_idx[i] += stride * blocksize;
          if (end_idx[i] > global_dimensions[i] - 1) {
            end_idx[i] = global_dimensions[i] - 1;
          }
        }
        block_interpolation(
            decData, block.get_global_index(), end_idx, PB_recover,
            interpolators[interpolator_id], direction_sequence_id, stride);
      }
    }
    quantizer.postdecompress_data();

    // postprocess to remove artifacts

    //            timer.stop("Interpolation Decompress");
    // std::cout << "counter = " << counter << std::endl;



    return decData;
  }

  // compress given the error bound
  // uchar *compress(const Config &conf, T *data, size_t &compressed_size)
  // {
  //   return compress(conf, data, compressed_size, true);
  // }

  uchar *compress(const Config &conf, T *data, size_t &compressed_size, bool tuning = false)
  {
//     std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
//     blocksize = conf.interpBlockSize;
//     interpolator_id = conf.interpAlgo;
//     direction_sequence_id = conf.interpDirection;
//     // assign additional variable



//     quant_pred_on = conf.quantization_prediction_on;
//     // if(tuning == true)
//     // {
//     //   quant_pred_on = false;
//     // }
//     quant_pred_start_level = conf.quantization_prediction_start_level;
//     // error smoothing



//     orig_data_ptr = (const T *) conf.PASS_DATA.original_data_prt;

//     init();
//     // original_variance = SZ::data_variance(data, num_elements); 
//     // std::cout << "original variance = " << original_variance << std::endl;
//     // For data range check.
//     auto orig_min_max = std::minmax_element(data, data + num_elements);
//     original_min = *orig_min_max.first;
//     original_max = *orig_min_max.second;
//     std::cout << "original max " << original_max << std::endl;
//     std::cout << "original min " << original_min << std::endl;
    
//     original_range = original_max - original_min;

//     Timer timer;

//     quant_inds.reserve(num_elements);
//     size_t interp_compressed_size = 0;


//     double eb_input = quantizer.get_eb();
//     double eb_final;  // eb for the highest level
//     // double eb_reduction_factor;
//     if (interpolator_id == 0) {
//       eb_final = eb_input /
//                  pow(linear_interp_eb_factor, (interpolation_level - 1) * N);
//       eb_reduction_factor = pow(linear_interp_eb_factor, N);
//     }
//     else {
//       eb_final = eb_input /
//                  pow(cubic_interp_eb_factor, (interpolation_level - 1) * N);
//       eb_reduction_factor = pow(cubic_interp_eb_factor, N);
//     }

//     // quant_inds.push_back(quantizer.quantize_and_over*data, 0));
//     quantize(0, *data, 0);

//     // Timer timer;
//     timer.start();

//     // double reduction_factor;
//     // double real_eb_ratio;
//     // if( interpolators[interpolator_id] == "linear")
//     // {
//     //     reduction_factor = sqrt(27/8);
//     // }
//     // else
//     // {
//     //     reduction_factor = sqrt(4.462681);
//     // }
//     // real_eb_ratio = pow(1/reduction_factor, interpolation_level-1);


//     for (uint level = interpolation_level;
//          level > 0 && level <= interpolation_level; level--) {
//       // direction_sequence_id = direction_choice(mt);
//       // std::cout << "direction_sequence_id = " << direction_sequence_id <<
//       // std::endl; for (int i = 0; i < N; i++)
//       // {
//       //     std::cout << dimension_sequences[direction_sequence_id][i] << " ";
//       // }
//       // std::cout << std::endl;
//       // if (level >= 3) {
//       //     quantizer.set_eb(eb_input * eb_ratio);
//       // } else {
//       //     quantizer.set_eb(eb_input);
//       // }
     
      
//       current_level = level;
//       current_base_eb = eb_final;
//       quantizer.set_eb(eb_final);
//       eb_final *= eb_reduction_factor;

//       // linear increase 
//       // current_level = level; 
//       // double current_eb = eb_input + (1e-32 - eb_input)/(interpolation_level -1)*(level-1);
//       // quantizer.set_eb(current_eb);

//       // exponent 
//       // current_level = level;
//       // double current_eb = eb_input * std::exp(1) * std::exp(-1.0*level);
//       // quantizer.set_eb(eb_input/(level*level*level));

//       // if(level >1)
//       // {
//       //   quantizer.set_eb(1e-32*eb_input);
//       // }
//       // else {
//       //   quantizer.set_eb(eb_input);
//       // }


//       // if (level <= 2) {
//       //   quantizer.set_eb(eb_input);
//       // }

//       size_t stride = 1U << (level - 1);

//       // printf("blocksize = %d\n", blocksize);
//       {
//       auto inter_block_range =
//           std::make_shared<SZ::multi_dimensional_range<T, N>>(
//               data, std::begin(global_dimensions), std::end(global_dimensions),
//               blocksize * stride, 0);

//       auto inter_begin = inter_block_range->begin();
//       auto inter_end = inter_block_range->end();

//       for (auto block = inter_begin; block != inter_end; ++block) {
//         auto end_idx = block.get_global_index();
//         for (int i = 0; i < N; i++) {
//           end_idx[i] += blocksize * stride;
//           if (end_idx[i] > global_dimensions[i] - 1) {
//             end_idx[i] = global_dimensions[i] - 1;
//           }
//         }
//         block_interpolation(
//             data, block.get_global_index(), end_idx, PB_predict_overwrite,
//             interpolators[interpolator_id], direction_sequence_id, stride);
//       }
//       }

//     }

//     assert(quant_inds.size() <= num_elements);

//     encoder.preprocess_encode(quant_inds, 0);
//     size_t bufferSize = 1.5 * (quantizer.size_est() + encoder.size_est() +
//                                sizeof(T) * quant_inds.size());

//     // TODO: change to smart pointer here
//     uchar *buffer = new uchar[bufferSize];
//     uchar *buffer_pos = buffer;

//     std::cout << "bufferSize = " << bufferSize << std::endl; 

//     write(global_dimensions.data(), N, buffer_pos);
//     write(blocksize, buffer_pos);
//     write(interpolator_id, buffer_pos);
//     write(direction_sequence_id, buffer_pos);

//     write(original_max, buffer_pos);
//     write(original_min, buffer_pos);


//     quantizer.save(buffer_pos);
//     quantizer.postcompress_data();

//     timer.start();
//     encoder.save(buffer_pos);
//     encoder.encode(quant_inds, buffer_pos);
//     encoder.postprocess_encode();
//     //            timer.stop("Coding");
//     assert(buffer_pos - buffer < bufferSize);
    
//     // writefile("compressed.dat",data, num_elements);
//     timer.start();
//     uchar *lossless_data =
//         lossless.compress(buffer, buffer_pos - buffer, compressed_size);

//     lossless.postcompress_data(buffer);




//       // free the error compression ptr 



//     // post process
//     // writefile("compressed.dat", data, num_elements);
//     // if(N==3)
//     // {
//     //   std::cout << "3D post process" << std::endl;
//     //   compensation_3d2(
//     //   data, aux_quant_inds.data(), conf.dims.data(), 0,quantizer.get_eb());
//     // }

//     // writefile("post_compressed.dat", data, num_elements);

//     // writefile("rand.dat", rand_collector.data(), num_elements);
//     // writefile("pred_noise.dat", my_pred_noise.data(), num_elements);
//     // writefile("error.dat", error_recorder.data(), num_elements);
//     // writefile("quant_inds.compress.dat", quant_inds.data(),
//     // quant_inds.size());

//     //   writefile("flushed_block_id.compress.dat", flushed_block_id.data(),
//     //           num_detection_block);
//     // writefile("significant_block_id.compress.dat",
//     //           significant_block_id.data(), num_detection_block);
//     // writefile("flushed_block.compress.dat", flushed_block.data(),
//     //           num_elements);
//     // writefile("significant_block.compress.dat", significant_block.data(),
//     // num_elements);


// #ifdef SZ_ERROR_ANALYSIS
//     writefile("pred.dat", my_pred.data(), num_elements);
//     writefile("quant.dat", aux_quant_inds.data(), num_elements);
//     writefile("quant_processed.dat", my_quant_inds_copy.data(), num_elements);
//     writefile("decompressed.dat", data, num_elements);
//     writefile("level.dat", my_level.data(), num_elements);
//     writefile(
//         "interp_direction.dat", my_interp_direction.data(), num_elements);
//     writefile(
//         "compensation_label.int32", my_compensation_label.data(),
//         my_compensation_label.size());
//     std::cout << "[ANALYSIS COMPILATION MODE]" << std::endl;

//     // Try huffman encoding on inorder quantization and level-wise quantization
//     // integers
//     // SZ::HuffmanEncoder<int> huff_coding = SZ::HuffmanEncoder<int>();

//     // uchar *test_buffer = new uchar[bufferSize];
//     // uchar *test_buffer_pos = test_buffer;
//     // huff_coding.preprocess_encode(aux_quant_inds, 0);
//     // huff_coding.save(test_buffer_pos);
//     // huff_coding.encode(aux_quant_inds, test_buffer_pos);
//     // huff_coding.postprocess_encode();
//     // size_t comsize = 0;
//     // uchar *lossless_data2 =
//     //     lossless.compress(test_buffer, test_buffer_pos - test_buffer,
//     //     comsize);
//     // // lossless.postcompress_data(test_buffer);

//     // std::cout << "[inorder]comsize = " << comsize << std::endl;
//     // free(test_buffer);

//     // // try it on level-wise quantization integers
//     // SZ::HuffmanEncoder<int> huff_coding2 = SZ::HuffmanEncoder<int>();

//     // uchar *test_buffer2 = new uchar[bufferSize];
//     // uchar *test_buffer_pos2 = test_buffer2;
//     // huff_coding2.preprocess_encode(quant_inds, 0);
//     // huff_coding2.save(test_buffer_pos2);
//     // huff_coding2.encode(quant_inds, test_buffer_pos2);
//     // huff_coding2.postprocess_encode();

//     // size_t comsize2 = 0;
//     // uchar *lossless_data3 = lossless.compress(
//     //     test_buffer2, test_buffer_pos2 - test_buffer2, comsize2);
//     // // lossless.postcompress_data(test_buffer2);

//     // std::cout << "[level-wise]comsize = " << comsize2 << std::endl;

//     // std::cout << "[level-wise] huffman encoding size = " << test_buffer_pos2
//     // - test_buffer2 << std::endl;
//     // free(test_buffer2);

//     // writefile("quant_level.dat", quant_inds.data(), quant_inds.size());

// #endif

//     compressed_size += interp_compressed_size;
    return nullptr;
    // return lossless_data;
  }

 private:
  enum PredictorBehavior { PB_predict_overwrite, PB_predict, PB_recover };

  void init()
  {
    assert(
        blocksize % 2 == 0 &&
        "Interpolation block size should be even numbers");
    num_elements = 1;
    interpolation_level = -1;
    for (int i = 0; i < N; i++) {
      if (interpolation_level < ceil(log2(global_dimensions[i]))) {
        interpolation_level = (uint)ceil(log2(global_dimensions[i]));
      }
      num_elements *= global_dimensions[i];
    }

    dimension_offsets[N - 1] = 1;
    for (int i = N - 2; i >= 0; i--) {
      dimension_offsets[i] =
          dimension_offsets[i + 1] * global_dimensions[i + 1];
      std::cout << "dimension_offsets[" << i << "] = " << dimension_offsets[i]
                << std::endl;
    }

    dimension_sequences = std::vector<std::array<int, N>>();
    auto sequence = std::array<int, N>();
    for (int i = 0; i < N; i++) { sequence[i] = i; }
    do {
      dimension_sequences.push_back(sequence);
    } while (std::next_permutation(sequence.begin(), sequence.end()));

    // post process utils 
    aux_quant_inds.resize(num_elements,0);

#ifdef SZ_ERROR_ANALYSIS
    my_level.resize(num_elements);
    my_quant_inds_copy.resize(num_elements);
    my_pred.resize(num_elements);
    my_pred[0] = 0;
    my_level[0] = interpolation_level;
    my_pred_noise.resize(num_elements, 0);
    my_interp_direction.resize(num_elements, 0);
    my_compensation_label.resize(num_elements, 0);
#endif

  }

  inline void range_check(T &d)
  {
    if (d > original_max) { d = original_max; }
    if (d < original_min) { d = original_min; }
  }

  inline int backward_compensate_pred(
      size_t idx, size_t offset1, size_t offset2, T &pred,
      const double compensation)
  {
    // return 0;
    // compensation = ? * eb; 0.5, 1.1, 1.5, 2
    //    0   1  2
    // 0 *A  *B *C
    // 1 *D  *E *F
    // 2 *G  *H  X
    // one layer is enough
    int A = idx - offset1 * 2 - offset2 * 2;
    int B = idx - offset1 * 2 - offset2;
    int C = idx - offset1 * 2;
    int D = idx - offset1 - offset2 * 2;
    int E = idx - offset1 - offset2;
    int F = idx - offset1;
    int G = idx - offset2 * 2;
    int H = idx - offset2;
    if (aux_quant_inds[E] == 0 || aux_quant_inds[F] == 0 ||
        aux_quant_inds[H] == 0) {
      return 0;
    }
    // int quant_A = aux_quant_inds[A] - quantizer.get_radius();
    // int quant_B = aux_quant_inds[B] - quantizer.get_radius();
    // int quant_C = aux_quant_inds[C] - quantizer.get_radius();
    // int quant_D = aux_quant_inds[D] - quantizer.get_radius(); // be cautious about out-of-bound access
    int quant_E = aux_quant_inds[E] - quantizer.get_radius();
    int quant_F = aux_quant_inds[F] - quantizer.get_radius();
    // int quant_G = aux_quant_inds[G] - quantizer.get_radius();
    int quant_H = aux_quant_inds[H] - quantizer.get_radius();

    int quant_compensate = 0;
    if (quant_H > 0 && quant_F > 0) {
      quant_compensate = (quant_H + quant_F - quant_E);
      pred += quant_compensate * compensation;
    }
    else if (quant_H < 0 && quant_F < 0) {
      quant_compensate = (quant_H + quant_F - quant_E);
      pred+= quant_compensate * compensation;
    }
    else {
      return 0;
    }
    return quant_compensate;
  }



  inline void quantize(size_t idx, T &d, T pred)
  {
    T d_copy = d;
       
      quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
    // range_check(d);
    aux_quant_inds[idx] = quant_inds.back();
#ifdef SZ_ERROR_ANALYSIS
    my_level[idx] = current_level;
    my_quant_inds_copy[idx] = quant_inds.back();
    // my_pred[idx] = pred;
    my_interp_direction[idx] = my_current_interp_direction;
#endif
  }

  inline void recover(size_t idx, T &d, T pred)
  {
    d = quantizer.recover(pred, quant_inds[quant_index++]);  
    // range_check(d);
    aux_quant_inds[idx] = quant_inds[quant_index-1];
  };

  double block_interpolation_1d(
      T *data, size_t begin, size_t end, size_t stride,
      const std::string &interp_func, const PredictorBehavior pb,
      bool quant_pred = false, size_t offset1 = 0, size_t offset2 = 0, bool is_left_boundary = true,
      bool use_fft_interp = false)
  {
    quant_pred = (quant_pred_on == true) && (is_left_boundary == false) && (current_level <= quant_pred_start_level);
    size_t n = (end - begin) / stride + 1;
    if (n <= 1) { return 0; }
    double predict_error = 0;

    int quant_compensation = 0;

    size_t stride3x = 3 * stride;
    size_t stride5x = 5 * stride;

    if (interp_func == "linear" || n < 5) {
      if (pb == PB_predict_overwrite) {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          T d_copy = *d;
          T pred = interp_linear(*(d - stride), *(d + stride));

          #ifdef SZ_ERROR_ANALYSIS
          my_pred[d - data] = pred;
          #endif

          if (quant_pred == true) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, compensation);            
          }

          quantize(d - data, *d, pred);
          if(quant_inds.back()==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          }
          #ifdef SZ_ERROR_ANALYSIS
          my_compensation_label[d - data] = quant_compensation;
          #endif

        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          T pred = interp_linear(*(d - stride), *(d - stride));
          #ifdef SZ_ERROR_ANALYSIS
          my_pred[d - data] = pred;
          #endif
          if ( quant_pred == true) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, compensation);
            // backward_compensate(d - stride-data, d - stride-data, pred,
            // error_recorder.data(), compensation, tol);
          }
          quantize(d - data, *d, pred);
          if(quant_inds.back()==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          }
          #ifdef SZ_ERROR_ANALYSIS
          my_compensation_label[d - data] = quant_compensation;
          #endif
          // quantize(d - data, *d, *(d - stride));

          // if (n < 4) {
          //     quantize(d - data, *d, *(d - stride));
          // } else {
          //     quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d -
          //     stride)));
          // }
        }
      }
      else {
        for (size_t i = 1; i + 1 < n; i += 2) {
          T *d = data + begin + i * stride;
          T pred = interp_linear(*(d - stride), *(d + stride));
          // pred prediction 
          if(quant_pred == true)
          {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation = region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
          }

          recover(d - data, *d, pred);

          if(quant_inds[quant_index-1]==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds[quant_index-1] + quant_compensation;
          }

        }
        if (n % 2 == 0) {
          T *d = data + begin + (n - 1) * stride;
          T pred = *(d - stride); 
          if(quant_pred == true)
          {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation = region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
          }
          recover(d - data, *d, pred);
          if(quant_inds[quant_index-1]==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds[quant_index-1] + quant_compensation;
          }
          // if (n < 4) {
          //     recover(d - data, *d, *(d - stride));
          // } else {
          //     recover(d - data, *d, interp_linear1(*(d - stride3x), *(d -
          //     stride)));
          // }
        }
      }
    }
    else {
      if (pb == PB_predict_overwrite) {
        double dafault_eb = quantizer.get_eb();
        T *d;
        size_t i;

        d = data + begin + stride;
        T pred = interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x));

        if(quant_pred == true)
        {
          double compensation = region_error_control_eb_compensation * quantizer.get_eb();
          quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
        }
        quantize(d - data, *d, pred);
        aux_quant_inds[d - data] = quant_inds.back();
        if(quant_inds.back()==0)
        {
          aux_quant_inds[d - data] = 0;
        }
        else {
          aux_quant_inds[d - data] = quant_inds.back() + quant_compensation;
        }



        #ifdef SZ_ERROR_ANALYSIS
        my_pred[d - data] = pred;
        #endif

        d = data + begin + stride;
        for (i = 3; i + 3 < n; i += 2) {
          d = data + begin + i * stride;
          T pred = interp_cubic(
              *(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x));
          #ifdef SZ_ERROR_ANALYSIS
          my_pred[d - data] = pred;
          #endif
          if ( quant_pred == true ) {
            double compensation = region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
          }
          quantize(d - data, *d, pred);
          if(quant_inds.back()==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          }
          #ifdef SZ_ERROR_ANALYSIS
          my_compensation_label[d - data] = quant_compensation;
          #endif
        }

        d = data + begin + i * stride;
        pred = interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride));
        if(quant_pred == true)
        {
          double compensation = region_error_control_eb_compensation * quantizer.get_eb();
          quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
        }
        quantize(d - data, *d, pred);
        if(quant_inds.back()==0)
        {
          aux_quant_inds[d - data] = 0;
        }
        else {
          aux_quant_inds[d - data] = quant_inds.back() + quant_compensation;
        }
        
        #ifdef SZ_ERROR_ANALYSIS
        my_pred[d - data] = pred;
        #endif

        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;
          // quantize(d - data, *d, *(d - stride));
          T pred = *(d - stride);

          #ifdef SZ_ERROR_ANALYSIS
          my_pred[d - data] = pred;
          #endif

          if ( quant_pred == true) {
            double tol = quantizer.get_eb() / linear_interp_eb_factor;
            double compensation =
                region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(
                d - data, offset1, offset2, pred, compensation);
            // backward_compensate(d - stride-data, d - stride-data, pred,
            // error_recorder.data(), compensation, tol);
          }
          quantize(d - data, *d, pred);
          if(quant_inds.back()==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds.back() + quant_compensation;
          }
          #ifdef SZ_ERROR_ANALYSIS
          my_compensation_label[d - data] = quant_compensation;
          #endif
          quantizer.set_eb(dafault_eb);
        }
      }
      else {
        T *d;
        size_t i;

        d = data + begin + stride;        
        recover(
            d - data, *d,
            interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x))); // beginning

        for (i = 3; i + 3 < n; i += 2) {
          d = data + begin + i * stride;
          T pred = interp_cubic(
              *(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x));
          if(quant_pred == true)
          {
            double compensation = region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
          }
          recover(d - data, *d, pred);
          if(quant_inds[quant_index-1]==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds[quant_index-1] + quant_compensation;
          }
        }
        d = data + begin + i * stride;
        recover(
            d - data, *d,
            interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride))); // end

        if (n % 2 == 0) {
          d = data + begin + (n - 1) * stride;
          T pred = *(d - stride);
          if(quant_pred == true)
          {
            double compensation = region_error_control_eb_compensation * quantizer.get_eb();
            quant_compensation = backward_compensate_pred(d - data, offset1, offset2, pred, compensation);
          }
          recover(d - data, *d, pred); // last element on the line sample
          if(quant_inds[quant_index-1]==0)
          {
            aux_quant_inds[d - data] = 0;
          }
          else {
            aux_quant_inds[d - data] = quant_inds[quant_index-1] + quant_compensation;
          }
        } 
      }
    }

    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 1, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    return block_interpolation_1d(
        data, begin[0], end[0], stride, interp_func, pb);
  }

  template <uint NN = N>
  typename std::enable_if<NN == 2, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    double predict_error = 0;
    size_t stride2x = stride * 2;
    // auto default_eb = quantizer.get_eb();

    // quantizer.set_eb(default_eb * c1);
    const std::array<int, N> dims = dimension_sequences[direction];
    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                            j * dimension_offsets[dims[1]];
      predict_error += block_interpolation_1d(
          data, begin_offset,
          begin_offset +
              (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
          stride * dimension_offsets[dims[0]], interp_func, pb);
    }
    // quantizer.set_eb(default_eb);
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      size_t begin_offset = i * dimension_offsets[dims[0]] +
                            begin[dims[1]] * dimension_offsets[dims[1]];
      predict_error += block_interpolation_1d(
          data, begin_offset,
          begin_offset +
              (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
          stride * dimension_offsets[dims[1]], interp_func, pb);
    }
    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 3, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    double predict_error = 0;
    size_t stride2x = stride * 2;

    auto default_eb = quantizer.get_eb();
    // quantizer.set_eb(default_eb * c2);

#ifdef SZ_ERROR_ANALYSIS
    my_current_interp_direction = 1;
#endif

    const std::array<int, N> dims = dimension_sequences[direction];
    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];
        // bool high_level_predict = (j%4==0) && (k%4==0);
        // high_level_predict = 0;
        bool is_left_boundary = (j==0) || (k==0);
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
            stride * dimension_offsets[dims[0]], interp_func, pb, quant_pred_on,
            dimension_offsets[dims[1]] * (stride2x),
            dimension_offsets[dims[2]] * (stride2x), is_left_boundary);
        quantizer.set_eb(default_eb);
      }
    }



    // quantizer.set_eb(default_eb * c1);
    // if(current_level ==1) quantizer.set_eb(default_eb*0.5);
#ifdef SZ_ERROR_ANALYSIS
    my_current_interp_direction = 2;
#endif
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              begin[dims[1]] * dimension_offsets[dims[1]] +
                              k * dimension_offsets[dims[2]];
        // if( i%stride2x ==0 || k%stride2x ==0)
        // {
        //   quantizer.set_eb(default_eb*0.8);
        // }
        bool is_left_boundary = (i==0) || (k==0); 
        // bool high_level_predict = (i%2==0) && (k%4==0);
        // high_level_predict = 0;
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
            stride * dimension_offsets[dims[1]], interp_func, pb, quant_pred_on,
            dimension_offsets[dims[0]] * (stride),
            dimension_offsets[dims[2]] * (stride2x),is_left_boundary);
        quantizer.set_eb(default_eb);
      }
    }
  

#ifdef SZ_ERROR_ANALYSIS
    my_current_interp_direction = 3;
#endif
    // if(current_level ==1) quantizer.set_eb(default_eb*0.5);
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      // ) std::cout <<"i = " << i << std::endl;
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        // if (i==0 && current_level==1) std::cout <<"j = " << j << std::endl;
        size_t begin_offset = i * dimension_offsets[dims[0]] +
                              j * dimension_offsets[dims[1]] +
                              begin[dims[2]] * dimension_offsets[dims[2]];
        // std::vector<size_t> cords =
        // interp_level_calculator.get_coordinates(begin_offset);
        // std::cout<<"cords = "<< cords[0] << " " << cords[1] << " " <<
        // cords[2] << std::endl;
        // if( i%stride2x ==1 && j%stride2x ==1)
        // {
        //   quantizer.set_eb(default_eb*0.65);
        // }
        bool is_left_boundary = (i==0) || (j==0); 
        bool use_fft_interp = false ; 
        predict_error += block_interpolation_1d(
            data, begin_offset,
            begin_offset +
                (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
            stride * dimension_offsets[dims[2]], interp_func, pb, quant_pred_on,
            dimension_offsets[dims[0]] * stride,
            dimension_offsets[dims[1]] * stride, is_left_boundary, use_fft_interp);
        quantizer.set_eb(default_eb);
        // if this is decompress and the level 1, we need to special treatment
      }
    }
    quantizer.set_eb(default_eb);
    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 4, double>::type block_interpolation(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    double predict_error = 0;
    size_t stride2x = stride * 2;
    max_error = 0;
    const std::array<int, N> dims = dimension_sequences[direction];
    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0);
         j <= end[dims[1]]; j += stride2x) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
             t <= end[dims[3]]; t += stride2x) {
          size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] +
                                j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
              stride * dimension_offsets[dims[0]], interp_func, pb);
        }
      }
    }
    max_error = 0;
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0);
           k <= end[dims[2]]; k += stride2x) {
        for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
             t <= end[dims[3]]; t += stride2x) {
          size_t begin_offset = i * dimension_offsets[dims[0]] +
                                begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
              stride * dimension_offsets[dims[1]], interp_func, pb);
        }
      }
    }
    max_error = 0;
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
             t <= end[dims[3]]; t += stride2x) {
          size_t begin_offset = i * dimension_offsets[dims[0]] +
                                j * dimension_offsets[dims[1]] +
                                begin[dims[2]] * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[2]] - begin[dims[2]]) * dimension_offsets[dims[2]],
              stride * dimension_offsets[dims[2]], interp_func, pb);
        }
      }
    }

    max_error = 0;
    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0);
         i <= end[dims[0]]; i += stride) {
      for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0);
           j <= end[dims[1]]; j += stride) {
        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0);
             k <= end[dims[2]]; k += stride) {
          size_t begin_offset = i * dimension_offsets[dims[0]] +
                                j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
          predict_error += block_interpolation_1d(
              data, begin_offset,
              begin_offset +
                  (end[dims[3]] - begin[dims[3]]) * dimension_offsets[dims[3]],
              stride * dimension_offsets[dims[3]], interp_func, pb);
        }
      }
    }
    return predict_error;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 1, double>::type block_post_process(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    return 0;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 2, double>::type block_post_process(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    return 0;
  }

  template <uint NN = N>
  typename std::enable_if<NN == 3, double>::type block_post_process(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    size_t stride2x = stride * 2;
    auto default_eb = quantizer.get_eb();
    const std::array<int, N> dims = dimension_sequences[direction];

      compensation_3d(data,aux_quant_inds.data(),
            begin.data(), end.data(), dims,
            dimension_offsets,
            0, stride, 
            1, 2, stride2x, stride2x, quantizer.get_eb(),
            quantizer.get_radius());
    

      compensation_3d(data,aux_quant_inds.data(),
            begin.data(), end.data(), dims,
            dimension_offsets,
            1,stride, 
            0, 2, stride, stride2x, quantizer.get_eb(),
            quantizer.get_radius());
    
  
    
      compensation_3d(data,aux_quant_inds.data(),
            begin.data(), end.data(), dims,
            dimension_offsets,
            2,stride, 
            0, 1, stride, stride, quantizer.get_eb(),
            quantizer.get_radius());

      // std::cout<<"quantizer radius = " << quantizer.get_radius() << std::endl;      
      writefile("first_compressed.dat", data, num_elements);
    
    return 0;
  }


  template <uint NN = N>
  typename std::enable_if<NN == 4, double>::type block_post_process(
      T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
      const PredictorBehavior pb, const std::string &interp_func,
      const int direction, size_t stride = 1)
  {
    return 0;
  }


  template <uint NN = N>
  typename std::enable_if<NN == 1, uchar*>::type error_sample_compress(const T* orig_data_ptr, T* data, 
        size_t& compressed_error_size, int direction_sequence_id,  
        const std::vector<std::array<int, 1>> &dimension_sequences, std::array<size_t,1>& global_dimensions, 
          std::array<size_t, 1>& dimension_offsets, uchar *&buffer_pos) 
  {
    return nullptr;

  }

  template <uint NN = N>
  typename std::enable_if<NN == 2, uchar*>::type error_sample_compress(const T* orig_data_ptr, T* data, 
        size_t& compressed_error_size, int direction_sequence_id,  
        const std::vector<std::array<int, 2>> &dimension_sequences, std::array<size_t,2>& global_dimensions, 
          std::array<size_t, 2>& dimension_offsets, uchar *&buffer_pos) 
  {
    return nullptr;
  }


  template <uint NN = N>
  typename std::enable_if<NN == 3, uchar*>::type error_sample_compress(const T* orig_data_ptr, T* data, 
        size_t& compressed_error_size, int direction_sequence_id,  
        const std::vector<std::array<int, 3>> &dimension_sequences, std::array<size_t,3>& global_dimensions, 
          std::array<size_t, 3>& dimension_offsets, uchar *&buffer_pos) 
  {


  }

  template <uint NN = N>
  typename std::enable_if<NN == 4, uchar* >::type error_sample_compress(const T* orig_data_ptr, T* data, 
        size_t& compressed_error_size, int direction_sequence_id,  
        const std::vector<std::array<int, 4>> &dimension_sequences, std::array<size_t,4>& global_dimensions, 
          std::array<size_t, 4>& dimension_offsets, uchar *&buffer_pos) 
  {
    return nullptr;
  }




  // at the lowest level, we construct a map to store the representitive's
  // quant index this is a bit map 1 for non-zero quant index, 0 otherwise only
  // for slices with

  int interpolation_level = -1;
  uint blocksize;
  int interpolator_id;
  double eb_ratio = 0.5;
  std::vector<std::string> interpolators = {"linear", "cubic"};
  std::vector<int> quant_inds;
  size_t quant_index = 0;  // for decompress
  double max_error;
  Quantizer quantizer;
  Encoder encoder;
  Lossless lossless;
  size_t num_elements;
  std::array<size_t, N> global_dimensions;
  std::array<size_t, N> dimension_offsets;
  std::vector<std::array<int, N>> dimension_sequences;
  int direction_sequence_id;
  // added for artifact mitigation
  int current_level = 0;
  // double c = sqrt(4.4159889); // deprecated
  // double c1 = 1.0 / sqrt(1.640625); //deprecated
  // double c2 = 1.0 / 1.640625; //deprecated
  // double c3 = 1.0 / sqrt(4.4159889); //deprecated

  double linear_interp_eb_factor = sqrt(1.5);
  double cubic_interp_eb_factor = 1.2808688457449497;


  T original_max;    // for data range check
  T original_min;    // for data range check
  T original_range;  // for data range check
  double current_base_eb;
  double eb_reduction_factor = 1.0; 

  // original data copy
  const T* orig_data_ptr;
  double original_variance;

  // post_process utils 

  std::vector<int> aux_quant_inds;
  double region_error_control_eb_compensation = 2.0; // for quant prediction 
  bool quant_pred_on = true;
  int quant_pred_start_level = 2;
   

// Analysis utils;
// This is for conditional compilation not comments
// #ifndef SZ_ERROR_ANALYSIS
// #define SZ_ERROR_ANALYSIS
#ifdef SZ_ERROR_ANALYSIS
  std::vector<int> my_level;
  std::vector<int> my_quant_inds_copy;
  std::vector<T> my_pred;
  // std::vector<T> my_pred_noise;
  std::vector<T> my_pred_noise;
  std::vector<int> my_interp_direction;
  int my_current_interp_direction = 0;
  std::vector<int> my_compensation_label;

#endif
};
};  // namespace SZ

#endif
