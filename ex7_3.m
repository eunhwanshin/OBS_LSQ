%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 7.3 in
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 143-147
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
% NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%

clear;

table6_2

% cofactor matrix of observations
Q = 0.01 * eye(6);

r2 = x(1,1)^2 + x(2,1)^2;
t1 = x(1,1) * y(1,1) + x(2,1) * y(2,1);
t2 = x(1,1) * y(2,1) - x(2,1) * y(1,1);

f = zeros(4,1);
j = 0;
for i=2:3,
  j = j+1;
  f(j) = t1*x(1,i) - t2*x(2,i) - y(1,i)*r2;
  j = j+1;
  f(j) = t1*x(2,i) + t2*x(1,i) - y(2,i)*r2; 
end
f = -f

A = zeros(4,6);
A(1,1) =  y(1,1)*x(1,2) - y(2,1)*x(2,2) - 2*y(1,2)*x(1,1);
A(1,2) =  y(2,1)*x(1,2) + y(1,1)*x(2,2) - 2*y(1,2)*x(2,1);
A(1,3) =  x(1,1)*y(1,1) + x(2,1)*y(2,1);
A(1,4) = -x(1,1)*y(2,1) + x(2,1)*y(1,1);

A(2,1) = y(1,1)*x(2,2) + y(2,1)*x(1,2) - 2*y(2,2)*x(1,1);
A(2,2) = y(2,1)*x(2,2) - y(1,1)*x(1,2) - 2*y(2,2)*x(2,1);
A(2,3) = x(1,1)*y(2,1) - x(2,1)*y(1,1);
A(2,4) = x(1,1)*y(1,1) + x(2,1)*y(2,1);

A(3,1) = y(1,1)*x(1,3) - y(2,1)*x(2,3) - 2*y(1,3)*x(1,1);
A(3,2) = y(2,1)*x(1,3) + y(1,1)*x(2,3) - 2*y(1,3)*x(2,1);
A(3,5) = A(1,3);
A(3,6) = A(1,4);

A(4,1) = y(1,1)*x(2,3) + y(2,1)*x(1,2) - 2*y(2,3)*x(1,1);
A(4,2) = y(2,1)*x(2,3) - y(1,1)*x(1,3) - 2*y(2,3)*x(2,1);
A(4,5) = A(2,3);
A(4,6) = A(2,4);

Qe = A * Q * A';
We = inv(Qe);
k = We *f
v = Q * A' * k

fprintf('Adjusted x\n')
x = x + [v(1), v(3), v(5); v(2), v(4), v(6)]

r2 = x(1,1)^2 + x(2,1)^2;
t1 = x(1,1)*y(1,1) + x(2,1)*y(2,1);
t2 = x(1,1)*y(2,1) - x(2,1)*y(1,1);

fprintf('Adjusted (a, b) = (%.3f, %.3f)\n', t1/r2, t2/r2)

% a posteriori cofactor matrix of the estimated observations
I_Qvv = eye(6) - Q * A' * We * A;
Qll = I_Qvv * I_Qvv * Q;

fprintf('Reference variance = %.6f\n', k'*f/4)
