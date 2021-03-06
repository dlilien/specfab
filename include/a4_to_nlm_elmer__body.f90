!*  (a_1    a_3         a_5)
!*  (a_3    a_2         a_4)
!*  (a_5    a_4   1-a_1-a_2)

!* 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
a4_1133 = a2(1) - a4(1) - a4(3)
a4_1233 = a2(3) - a4(7) - a4(9)
a4_1333 = a2(5) - a4(6) - a4(5)
a4_2233 = a2(2) - a4(2) - a4(3)
a4_2333 = a2(4) - a4(4) - a4(8)
a4_3333 = 1.0 - 2.0 * a2(1) - 2.0 * a2(2) + 2.0 * a4(3) + a4(1) + a4(2)

nlm(I_l0)=((15.0) + (-70.0) + (63.0)*a4(1) + (126.0)*a4(3) + (126.0)*a4_1133 + (63.0)*a4(2) + (126.0)*a4_2233 + (63.0)*a4_3333)/(16.*Sqrt(Pi));
nlm(I_l2+0)=(Sqrt(7.5)*((-7.0)*a2(1) - (0,14)*a2(3) + (7.0)*a2(2) + (9.0)*a4(1) + (0,18)*a4(7) + (9.0)*a4_1133 + (0,18)*a4(9) + (0,18)*a4_1233 + (-9.0)*a4(2) + (-9.0)*a4_2233)*Sqrt(Pi**(-1.0)))/8.;
nlm(I_l2+1)=(Sqrt(7.5)*((-7.0)*a2(5) - (0,7)*a2(4) + (9.0)*(a4(6) + (0,1)*a4(4) + a4(5) + a4_1333 + (0,1)*a4(8) + (0,1)*a4_2333))*Sqrt(Pi**(-1.0)))/4.;
nlm(I_l2+2)=(((7.0)*a2(1) + (7.0)*a2(2) + (-14.0)*(1.0 - a2(1) - a2(2)) + (-9.0)*a4(1) + (-18.0)*a4(3) + (9.0)*a4_1133 + (-9.0)*a4(2) + (9.0)*a4_2233 + (18.0)*a4_3333)*Sqrt((5.0)*Pi**(-1.0)))/8.;
nlm(I_l2+3)=(Sqrt(7.5)*((7.0)*a2(5) - (0,7)*a2(4) + (-9.0)*a4(6) + (0,9)*(a4(4) + (0,1)*a4(5) + (0,1)*a4_1333 + a4(8) + a4_2333))*Sqrt(Pi**(-1.0)))/4.;
nlm(I_l2+4)=-0.125*(Sqrt(7.5)*((7.0)*a2(1) - (0,14)*a2(3) + (-7.0)*a2(2) + (-9.0)*a4(1) + (0,18)*a4(7) + (-9.0)*a4_1133 + (0,18)*a4(9) + (0,18)*a4_1233 + (9.0)*a4(2) + (9.0)*a4_2233)*Sqrt(Pi**(-1.0)));
nlm(I_l4+0)=(3*Sqrt(17.5)*(a4(1) + (0,4)*a4(7) + (-6.0)*a4(3) - (0,4)*a4(9) + a4(2))*Sqrt(Pi**(-1.0)))/16.;
nlm(I_l4+1)=(3*(a4(6) + (0,3)*a4(4) + (-3.0)*a4(5) - (0,1)*a4(8))*Sqrt((35.0)*Pi**(-1.0)))/8.;
nlm(I_l4+2)=(-3*Sqrt(2.5)*(a4(1) + (0,2)*a4(7) + (-6.0)*a4_1133 + (0,2)*a4(9) - (0,12)*a4_1233 + (-1.0)*a4(2) + (6.0)*a4_2233)*Sqrt(Pi**(-1.0)))/8.;
nlm(I_l4+3)=(-3*((3.0)*a4(6) + (0,3)*a4(4) + (3.0)*a4(5) + (-4.0)*a4_1333 + (0,3)*a4(8) - (0,4)*a4_2333)*Sqrt((5.0)*Pi**(-1.0)))/8.;
nlm(I_l4+4)=(3*((3.0)*a4(1) + (6.0)*a4(3) + (-24.0)*a4_1133 + (3.0)*a4(2) + (-24.0)*a4_2233 + (8.0)*a4_3333))/(16.*Sqrt(Pi));
nlm(I_l4+5)=(3*((3.0)*a4(6) - (0,3)*a4(4) + (3.0)*a4(5) + (-4.0)*a4_1333 - (0,3)*a4(8) + (0,4)*a4_2333)*Sqrt((5.0)*Pi**(-1.0)))/8.;
nlm(I_l4+6)=(-3*Sqrt(2.5)*(a4(1) - (0,2)*a4(7) + (-6.0)*a4_1133 - (0,2)*a4(9) + (0,12)*a4_1233 + (-1.0)*a4(2) + (6.0)*a4_2233)*Sqrt(Pi**(-1.0)))/8.;
nlm(I_l4+7)=(-3*(a4(6) - (0,3)*a4(4) + (-3.0)*a4(5) + (0,1)*a4(8))*Sqrt((35.0)*Pi**(-1.0)))/8.;
nlm(I_l4+8)=(3*Sqrt(17.5)*(a4(1) - (0,4)*a4(7) + (-6.0)*a4(3) + (0,4)*a4(9) + a4(2))*Sqrt(Pi**(-1.0)))/16.;
