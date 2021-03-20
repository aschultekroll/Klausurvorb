#include <complex>
#include <fstream>
#include <vector>
#include <mpi.h>
#include<iostream>
#include<chrono>
using std::chrono::system_clock;
using std::chrono::duration;

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

// writes an image in PBM format
void write_pbm(std::ostream &os, std::vector<float> image, int width,
               int height) {
  os << "P1\n";
  os << width << ' ' << height << '\n';
  for (auto v : image)
    os << v << ' ';
}

// generate a nice picture
int main(int argc, char* argv[]) {//jeweiliger Rank
  //int argc;
  //char** argv;
  // image dimensions
  int width = 1920;
  int height = 1080;
  // the result vector
  std::vector<float> image;
  float* rbuf;
  rbuf = (float *)malloc(width * height*sizeof(float));
  // main computation
  MPI_Init(&argc,&argv);
  MPI_Status stat;
  int num;
  int rank;    
  MPI_Comm_size(MPI_COMM_WORLD, &num); // Anzahl der Prozesse
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout<<"rank "<<rank<<" of "<<num<<std::endl;
  int newheight=height/num;
  //std::vector<float>mpiimage(width*height);
  float* sendbuf;
  sendbuf = (float *)malloc(width * newheight*sizeof(float));
  //int start=rank*newheight;
  //int end=(rank+1)*newheight;
  
  //Zeitpunkt vor der Berechnung
  auto start = system_clock::now();
  MPI_Scatter(rbuf,width * height , MPI_FLOAT, rbuf, width * newheight, MPI_FLOAT,0,MPI_COMM_WORLD);
  std::cout<< "berechnung abgeschlossen"<<std::endl;
  for (int i = 0; i < newheight; ++i)
    for (int j = 0; j < width; ++j)
      sendbuf[(i) * width + j] = evaluate<100>(
          {3 * (j / float(width)) - 2, 2 * ((i+(rank*newheight)) / float(height) - 0.5f)});
    
    std::cout<< "berechnung abgeschlossen"<<std::endl;
  //Ergebeniss in end Vektor schreiben
  MPI_Gather(sendbuf,width * newheight , MPI_FLOAT, rbuf, width * newheight, MPI_FLOAT,0,MPI_COMM_WORLD);
  auto end= system_clock::now();
  const double elapsed_seconds = duration<double>(end-start).count();
  /*for(int r=0;r<num;++r)
    if(rank==r){
      int start=rank*newheight;
      int end=(rank+1)*newheight;  
      for (int i = start; i < end; ++i){
        for (int j = 0; j < width; ++j)
            image[i * width + j] = mpiimage[i * width + j];
      }
    }*/
  // writing the result*/
  if(rank==0){
    for(int i=0;i<width*height;++i)
      image.push_back(rbuf[i]);
  std::ofstream ofs("mpi.pbm");
  write_pbm(ofs, image, width, height);
  std::cout<<elapsed_seconds<<std::endl;
  }
  MPI_Finalize();  
  
  return 0;
}