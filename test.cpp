#include "datagen.h"
#include "p4p_solver.h"
#include "p35p_solver.h"
#include "runFunctions.h"
#include "gtest/gtest.h"


// Then we define a test fixture class template.
template <class T>
class SumTest : public testing::Test {
protected:
    // The ctor calls the factory function to create a prime table
    // implemented by T.
    SumTest() {}

    ~SumTest() override {}

    // Note that we test an implementation via the base interface
    // instead of the actual implementation class.  This is important
    // for keeping the tests close to the real world scenario, where the
    // implementation is invoked via the base interface.  It avoids
    // got-yas where the implementation class has a method that shadows
    // a method with the same name (but slightly different argument
    // types) in the base interface, for example.
};

#if GTEST_HAS_TYPED_TEST

using testing::Types;

// Google Test offers two ways for reusing tests for different types.
// The first is called "typed tests".  You should use it if you
// already know *all* the types you are gonna exercise when you write
// the tests.

// To write a typed test case, first use
//
//   TYPED_TE
//   ST_SUITE(TestCaseName, TypeList);
//
// to declare it and specify the type parameters.  As with TEST_F,
// TestCaseName must match the test fixture name.

using MyTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(SumTest, MyTypes);
// The list of types we want to test.
// Then use TYPED_TEST(TestCaseName, TestName) to define a typed test,
// similar to TEST_F..

TYPED_TEST(SumTest, eqs) {
    ASSERT_EQ(5, 5);
}

// That's it!  Google Test will repeat each TYPED_TEST for each type
// in the type list specified in TYPED_TEST_SUITE.  Sit back and be
// happy that you don't have to define them multiple times.

#endif  // GTEST_HAS_TYPED_TEST