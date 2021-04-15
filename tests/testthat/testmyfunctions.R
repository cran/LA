
test_that("Generate a random OofA matrix with dimension n by k", {
  expect_equal(dim(rOofA(n=6,k=3)),c(6,3))
  expect_equal(dim(rOofA(n=9,k=4)),c(9,4))
})

test_that("Generate the pairwise-order (PWO) of an order-of-addition design", {
  expect_equal(dim(PWO(rOofA(k=3))),c(6,3))
  expect_equal(dim(PWO(rOofA(k=4))),c(24,6))
})

test_that("Generate the moment matrix of an order-of-addition design", {
  expect_equal(dim(MOM(rOofA(k=3))),c(4,4))
  expect_equal(dim(MOM(rOofA(k=4))),c(7,7))
})

test_that("Generate an optimal order-of-addition design with dimension n by k via LA", {
  expect_equal(dim(LA_OofA(n=6,k=3)),c(6,3))
  expect_equal(dim(LA_OofA(n=9,k=4)),c(9,4))
})

test_that("Generate an optimal LHD with dimension n by k via LA", {
  expect_equal(dim(LA_LHD(n=6,k=3)),c(6,3))
  expect_equal(dim(LA_LHD(n=9,k=4)),c(9,4))
})

test_that("Generate an optimal solution for function optimization via LA", {
  SDP=function(x){i=1:length(x);y=sum(abs(x)^(i=1));return(y)}
  CiT=function(x){x1=x[1];x2=x[2];
  y=-0.0001*(abs(sin(x1)*sin(x2)*exp(abs(100-sqrt(x1^2+x2^2)/pi)))+1)^0.1;return(y)}

  expect_equal(dim(LA_opt(of=SDP,lb=rep(-1,20),ub=rep(1,20),N=10,type="mini")),NULL)
  expect_equal(dim(LA_opt(of=CiT,lb=rep(-10,2),ub=rep(10,2),N=10,type="mini")),NULL)
})
