// ClusterShapeHitFilter test
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

#include <iostream>
#include <cassert>

namespace test {
namespace ClusterShapeHitFilterTest {
int test_limits()
{
    static const std::string pixelShapeFile("RecoPixelVertexing/PixelLowPtUtilities/data/pixelShape.par");
    static const std::string stripShapeFile("RecoPixelVertexing/PixelLowPtUtilities/data/stripShape.par");
    const cluster_shape::PixelLimitsCollection pixelLimits(pixelShapeFile);
    const cluster_shape::StripLimitsCollection stripLimits(stripShapeFile);

    static constexpr float out = 10e12;
    static constexpr float eps = 0.01;
    std::cout << "dump strip limits" << std::endl;
    for(int i = 0; i != cluster_shape::StripKeys::N + 1; ++i) {
        assert(!stripLimits[i].isInside(out));
        assert(!stripLimits[i].isInside(-out));
        std::cout << i << ": ";
        const float* p = stripLimits[i].data[0];
        if(p[1] < 1.e9) {
            assert(stripLimits[i].isInside(p[0] + eps));
            assert(stripLimits[i].isInside(p[3] - eps));
        }
        for(int j = 0; j != 4; ++j)
            std::cout << p[j] << ", ";
        std::cout << std::endl;
    }

    const std::pair<float, float> out1(out, out), out2(-out, -out);
    std::cout << "\ndump pixel limits" << std::endl;
    for(int i = 0; i != cluster_shape::PixelKeys::N + 1; ++i) {
        assert(!pixelLimits[i].isInside(out1));
        assert(!pixelLimits[i].isInside(out2));
        std::cout << i << ": ";
        const float* p = pixelLimits[i].data[0][0];
        if (p[1]<1.e9) {
            assert(pixelLimits[i].isInside(std::pair<float, float>(p[0] + eps, p[3] - eps)));
            assert(pixelLimits[i].isInside(std::pair<float, float>(p[5] - eps, p[6] + eps)));
        }
        for(int j = 0; j != 8; ++j)
            std::cout << p[j] << ", ";
        std::cout << std::endl;
    }

    return 0;
}
} // namespace ClusterShapeHitFilterTest
} // namespace test

int main()
{
    return test::ClusterShapeHitFilterTest::test_limits();
}
