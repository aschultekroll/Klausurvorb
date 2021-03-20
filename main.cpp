#include <complex> 
#include <fstream> 
#include <vector>

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
void write_pbm(std::ostream &os, std::vector<float> image, int width, int height) {
    os << "P1\n";
    os << width << ' ' << height << '\n'; for (auto v : image)
    os << v << ' '; }
    // generate a nice picture 
int main() {
    // image dimensions
    int width = 1920;
    int height = 1080;
    // the result vector
    std::vector<float> image(width * height); 
    // main computation
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j){
        image[i * width + j] = evaluate<100>(
            {3 * (j / float(width)) - 2, 2 * (i / float(height) - 0.5f)});
        }}
    // writing the result
    std::ofstream ofs("image.pbm"); 
    write_pbm(ofs, image, width, height); return 0;
}

