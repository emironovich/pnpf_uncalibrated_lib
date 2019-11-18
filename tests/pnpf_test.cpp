#include "test_helper.hpp"
#include "gtest/gtest.h"

using namespace pnpf;

template <typename Solver> struct ExpectedResults {
  static constexpr double solution_rate = 0.99;
  static constexpr double success_rate = 0.99;
};

template <typename Solver> class PnPTest : public testing::Test {
protected:
  PnPTest() {
    it_num = 1e+4; // number of iterations
    Solver solver;
    res = runFunction(solver, it_num);
  };
  int it_num;
  TestResult res;
  using Expected = ExpectedResults<Solver>;
  const double solution_rate_exp = Expected::solution_rate;
  const double success_rate_exp = Expected::success_rate;
};

using testing::Types;

using SolverTypes = ::testing::Types<P35PSolver<double>, P35PSolver<float>,
                                     P4PSolver<double>, P4PSolver<float>>;

template <> struct ExpectedResults<P35PSolver<float>> {
  static constexpr double solution_rate = 0.95;
  static constexpr double success_rate = 0.9;
};
template <> struct ExpectedResults<P4PSolver<float>> {
  static constexpr double solution_rate = 0.85;
  static constexpr double success_rate = 0.75;
};

TYPED_TEST_SUITE(PnPTest, SolverTypes);

TYPED_TEST(PnPTest, PnP) {
  EXPECT_GT((double)this->res.existSolutions / this->it_num,
            this->solution_rate_exp);
  EXPECT_GT((double)this->res.belowThreshold / this->it_num,
            this->success_rate_exp);
}
