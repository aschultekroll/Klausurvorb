#include <complex> 
#include <fstream> 
#include <vector>
#include <chrono>

#include <vector>
#include <boost/align/aligned_allocator.hpp> 
#include <omp.h>
#include<iostream>
// the number type we use here 
using ct = std::complex<float>;

        
// evaluates the main condition on an input
template <int MAX_ITER> float evaluate(const ct &p) {
    constexpr int LIMIT = 10;
    ct l(0, 0);
    for (int i = 0; i < MAX_ITER; ++i)
    l = l * l + p;
    return (std::norm(l) > LIMIT);
}



constexpr int ALIGNMENT = 32; 
template <typename T> 
using my_alloc=boost :: alignment :: aligned_allocator<T, ALIGNMENT>;

template <typename T> 
using aligned_vector = std :: vector<T, my_alloc<T>> ;



// writes an image in PBM format
void write_pbm(std::ostream &os, aligned_vector<float> image, int width, int height) {
    os << "P1\n";
    os << width << ' ' << height << '\n'; for (auto v : image)
    os << v << ' '; }
    // generate a nice picture 
int main() {
    // image dimensions
    int width = 1920;
    int height = 1080;
    // the result vector
    
    
    //...
    // allocate and initialize aligned storage
        //aligned_vector<float> a(N, 1.1f);
        //aligned_vector<float> b(N, 2.2f);
    // compute vectorized
    //float∗ ap = a.data(); 
    //float∗ bp = b.data();
    
    // float sum = 0.0f;
    // #pragma omp simd aligned(ap,bp:ALIGNMENT) reduction(+:sum) 
    // for (int i = 0; i < N; ++i) {
    //     sum += a [ i ] * b [ i ] ;
    // }
    
    aligned_vector<float> image(width * height); 
    float* imagep=image.data();
    // main computation
    
    auto start = std::chrono::system_clock::now();
    
    for (int i = 0; i < height; ++i) {
        #pragma omp simd aligned(imagep:ALIGNMENT)
        for (int j = 0; j < width; ++j){
        imagep[i * width + j] = evaluate<100>(
            {3 * (j / float(width)) - 2, 2 * (i / float(height) - 0.5f)});
        }}
    auto end = std::chrono::system_clock::now();
    const double elapsed_seconds = std::chrono::duration<double>(end - start).count();
    std::cout<<"Time: "<< elapsed_seconds<<std::endl;
    // writing the result
    std::ofstream ofs("image.pbm"); 
    write_pbm(ofs, image, width, height); return 0;
}

