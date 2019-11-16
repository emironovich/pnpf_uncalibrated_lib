#include "datagen.h"
#include "runFunctions.h"
#include "gtest/gtest.h"

template <class T> class PnPTest : public testing::Test {
protected:
  PnPTest() {
    it_num = 1e+4; // number of iterations
    P35PSolver<T> p35pSolver;
    P4PSolver<T> p4pSolver;
    p35pRes = runFunction(p35pSolver, it_num);
    p4pRes = runFunction(p4pSolver, it_num);
  };
  int it_num;
  TestResult p35pRes, p4pRes;
};

#if GTEST_HAS_TYPED_TEST

using testing::Types;

using PrecisionTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(PnPTest, PrecisionTypes);

TYPED_TEST(PnPTest, P35P) {
  ASSERT_GT((double)this->p35pRes.existSolutions / this->it_num,
            0.8); // todo: figure out constants
  ASSERT_GT((double)this->p35pRes.belowThreshold / this->it_num, 0.8);
}

TYPED_TEST(PnPTest, P4P) {
  ASSERT_GT((double)this->p4pRes.existSolutions / this->it_num, 0.8);
  ASSERT_GT((double)this->p4pRes.belowThreshold / this->it_num, 0.77);
}

#endif // GTEST_HAS_TYPED_TEST