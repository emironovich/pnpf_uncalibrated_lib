#include "datagen.h"
#include "runFunctions.h"
#include "gtest/gtest.h"
#include <utility>

template <class SolverClass> class PnPTest : public testing::Test {
protected:
  PnPTest() {
    it_num = 1e+4; // number of iterations
    SolverClass solver;
    res = runFunction(solver, it_num);
  };
  int it_num;
  TestResult res;
};

#if GTEST_HAS_TYPED_TEST

using testing::Types;

using PrecisionTypes = ::testing::Types<P35PSolver<double>, P35PSolver<float>,
                                        P4PSolver<double>, P4PSolver<float>>;
TYPED_TEST_SUITE(PnPTest, PrecisionTypes);

TYPED_TEST(PnPTest, PnP) {
  ASSERT_GT((double)this->res.existSolutions / this->it_num,
            0.8); // todo: figure out constants
  ASSERT_GT((double)this->res.belowThreshold / this->it_num, 0.8);
}

#endif // GTEST_HAS_TYPED_TEST