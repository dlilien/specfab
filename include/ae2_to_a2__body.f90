!*  (a_1    a_3         a_5)
!*  (a_3    a_2         a_4)
!*  (a_5    a_4   1-a_1-a_2)

a2(1, 1) = ae2(1)
a2(2, 2) = ae2(2)
a2(3, 3) = 1.0 - ae2(1) - ae2(2)
a2(1, 2) = ae2(3)
a2(2, 1) = ae2(3)
a2(3, 2) = ae2(4)
a2(2, 3) = ae2(4)
a2(3, 1) = ae2(5)
a2(1, 3) = ae2(5)
