test_that("Hard ranking loss is implemented correctly",{
y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
expect_equal(Rank()@risk(y,yhat),8/45)
expect_error(Rank()@risk(c(1,2,3,4),c(1,2,3)))
expect_equal(Rank()@risk(c(1,6,5,3,7,8,2),c(6,2,9,1,-2,4,5)),sum(sign(outer(c(1,6,5,3,7,8,2),c(1,6,5,3,7,8,2),function(x,z)z-x))-sign(outer(c(6,2,9,1,-2,4,5),c(6,2,9,1,-2,4,5),function(x,z)z-x))!=0)/42)
})