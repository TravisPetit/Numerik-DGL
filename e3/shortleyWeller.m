%% shortleyWeller
% Input
% C see @makeGrid
% H see @makeGrid
function [A, Ab] = shortleyWeller(C, H)
nd = max(C, [], 'all');
ai = [];
aj = [];
av = [];

for j = 1:nd
    ai = [ai ; j ; j      ; j      ; j      ; j];
    aj = [aj ; j ; C(1,j) ; C(2,j) ; C(3,j) ; C(4,j)];
    
    t0 = 2/(H(1,j) * H(2,j)) + 2/(H(3,j) * H(4,j));
    t1 = - H(1,j) * (H(1,j) + H(2,j));
    t2 = - H(2,j) * (H(1,j) + H(2,j));
    t3 = - H(3,j) * (H(3,j) + H(4,j));
    t4 = - H(4,j) * (H(3,j) + H(4,j));
    
    av = [av; t0; 2/t1; 2/t2; 2/t3; 2/t4];
end

ai_hat = ai(aj>=0);
aj_hat = aj(aj>=0);
av_hat = av(aj>=0);

ai_tilde = ai(aj<=0);
aj_tilde = aj(aj<=0);
av_tilde = av(aj<=0);

% size(ai)
% size(aj)
% size(av)
%
% size(ai_hat)
% size(aj_hat)
% size(av_hat)
% 
% size(ai_tilde)
% size(aj_tilde)
% size(av_tilde)
% 
% ai_hat
% aj_hat
% av_hat
% 
% ai_tilde
% aj_tilde
% av_tilde

A = sparse(ai_hat, aj_hat, av_hat);
Ab = sparse(ai_tilde, -aj_tilde, av_tilde);
end
