#include "datagen.h"
#include "p4p_solver.h"
#include "p35p_solver.h"
#include "runFunctions.h"
#include "gtest/gtest.h"

template <class T>
class PnPTest : public testing::Test {
protected:

    PnPTest() {
        it_num = 1e+3; //number of iterations
        p35pRes = runFunction("p35p", (T)1e-4, it_num); //todo: figure out constants
        p4pRes = runFunction("p4p", (T)1e-4, it_num);
        };
    int it_num;
    TestResult p35pRes, p4pRes;
};

#if GTEST_HAS_TYPED_TEST

using testing::Types;

using PrecisionTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(PnPTest, PrecisionTypes);

TYPED_TEST(PnPTest, P35P) {
    ASSERT_LT(this->p35pRes.existSolutions / this->it_num, 0.9); //todo: figure out constants
    ASSERT_LT(this->p35pRes.belowThreshold / this->it_num, 0.9);
}

TYPED_TEST(PnPTest, P4P) {
    ASSERT_LT(this->p4pRes.existSolutions / this->it_num, 0.9);
    ASSERT_LT(this->p4pRes.belowThreshold / this->it_num, 0.9);
}

#endif  // GTEST_HAS_TYPED_TEST