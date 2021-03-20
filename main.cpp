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

        
#include <mpi.h>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;







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
void write_pbm(std::ostream &os, std::vector<float> image, int width, int height) {
    os << "P1\n";
    os << width << ' ' << height << '\n'; for (auto v : image)
    os << v << ' '; }
    // generate a nice picture 
int main(int argc, char* argv[]) {
    int size, rank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status stat;

    std::cout<< "size: "<<size<<std::endl;
    std::cout<<"rank "<<rank<<" of "<<size<<std::endl;

    // image dimensions
    int width = 1920;
    int height = 1080;
    int newheight = 1000/size;
    // the result vector
    
    std::vector<float> image(width * height); 
    //float* imagep=image.data();
    float* rbuf;
    rbuf = (float *)malloc(width * height*sizeof(float));
    //std::vector<float> send_buf(newheight*width);
    float* sendbuf;
    sendbuf = (float *)malloc(width * newheight*sizeof(float));
    //std::vector<float> recv_buf(newheight*width);
    // main computation

    auto start = std::chrono::system_clock::now();
    MPI_Scatter(rbuf,width*height,MPI_FLOAT,rbuf,width*newheight,MPI_FLOAT,0,MPI_COMM_WORLD);
    for (int i = 0; i < newheight; ++i) {
        //#pragma omp simd aligned(imagep:ALIGNMENT)
        for (int j = 0; j < width; ++j){
        sendbuf[i * width + j] = evaluate<100>(
            {3 * (j / float(width)) - 2, 2 * ((i+(rank*newheight)) / float(height) - 0.5f)});
        }}
    auto end = std::chrono::system_clock::now();
    MPI_Gather(sendbuf,width*newheight,MPI_FLOAT,rbuf,width*newheight,MPI_FLOAT,0,MPI_COMM_WORLD);
    const double elapsed_seconds = std::chrono::duration<double>(end - start).count();
    std::cout<<"Time: "<< elapsed_seconds<<std::endl;

    if(rank==0){
    for(int i=0;i<width*height;++i){
      image.push_back(rbuf[i]);
    }
    // writing the result
    std::ofstream ofs("image.pbm"); 
    write_pbm(ofs, image, width, height); 
    MPI_Finalize();
    }
    return 0;
}


