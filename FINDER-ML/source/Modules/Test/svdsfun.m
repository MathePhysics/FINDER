function y = svdsfun(x,tflag,M,cx);

%disp(tflag);
%C = (cx' * cx) / M;
%y = C *x;

 y = (cx * x) ;
 y = cx' * y;
 y = y / M;
 
