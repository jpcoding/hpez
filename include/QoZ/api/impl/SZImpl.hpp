#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

#include "QoZ/def.hpp"
#include "QoZ/api/impl/SZDispatcher.hpp"
#include "QoZ/api/impl/SZImplOMP.hpp"
#include <cmath>
#include <memory>


template<class T, QoZ::uint N>
char *SZ_compress_impl(QoZ::Config &conf, const T *data, size_t &outSize) {
#ifndef _OPENMP
    conf.openmp=false;
#endif
    if (conf.openmp) {
        return SZ_compress_OMP<T, N>(conf, data, outSize);
    } else {
        // std::vector<T> dataCopy(data, data + conf.num);
        auto data_copy = std::make_shared<std::vector<T>>();
        data_copy->assign(data, data + conf.num);
        //std::cout<<"implstart"<<std::endl;
        auto output=SZ_compress_dispatcher<T, N>(conf, data_copy->data(), outSize);
        //std::cout<<"implend"<<std::endl;
        conf.PASS_DATA.processed_data_prt = std::static_pointer_cast<void>(data_copy);  
        return output;
    }
}


template<class T, QoZ::uint N>
void SZ_decompress_impl(QoZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
#ifndef _OPENMP
    conf.openmp=false;
#endif
    //std::cout<<"impl"<<conf.absErrorBound<<std::endl;
    if (conf.openmp) {
        SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
    }
}

#endif