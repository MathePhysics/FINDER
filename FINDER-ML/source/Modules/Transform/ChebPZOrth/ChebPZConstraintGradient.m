function [gradcon] = PZConstraintGradient(X, sub2, Datas, parameters, methods)


% First gradient constraint is just negative gradient of theta(1)
% gradcon1 = - sub2.dtheta_dy(1,:) * sub2.dy_dk;
% 
% % Second gradient constraint is just that of the Type II error bound
% gradcon2 = sub2.dfII_dgamma * sub2.dgamma_dtheta' * sub2.dtheta_dy * sub2.dy_dk;
% 
% gradcon = [gradcon1(:)' ; gradcon2(:)' ];

gradcon = (sub2.dg_dtheta * sub2.dtheta_dy * sub2.dy_dk)';

end