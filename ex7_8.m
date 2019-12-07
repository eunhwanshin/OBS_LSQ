%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 7.8 which solves 
% Example 6.2 using the technique of indirect observations.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 165-167.
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

% generates 
table6_2;

W = eye(6);

fprintf('=== Nonlinear Approach ===\n');

% approximate values
a = 1;
b = 2;

iter = 0;
corr = 1.0e+6;

while corr > 1.0e-6,

  iter = iter + 1;

  a2 = a^2;
  b2 = b^2;
  r2 = a2 + b2;
  iA = [a, b; -b, a] ./ r2;

  % Construct linearized model: v + B delta = f.

  f = zeros(6,1);
  B = zeros(6,2);
  B1 = [a2-b2, 2*a*b; -2*a*b, a2-b2] ./ r2^2;
  B2 = [2*a*b, b2-a2;  a2-b2, 2*a*b] ./ r2^2;
  for i=1:3,
    k1 = 2*i - 1;
    k2 = k1 + 1;
    f(k1:k2,1) = iA*y(:,i) - x(:,i);
    B(k1:k2,1) = B1 * y(:,i);
    B(k1:k2,2) = B2 * y(:,i);
  end

  % solve
  N = B'*W*B;
  iN = inv(N);
  t = B'*W*f;

  delta = iN * t;

  % update the estimate
  a = a + delta(1);
  b = b + delta(2);

  corr = norm(delta);

  fprintf('Iteration %d: a = %.6f, b = %.6f\n', iter, a, b);

end % while 

fprintf('=== Linear Approach ===\n');
fprintf('Let c = a/(a^2+b^2) and d = b/(a^2+b^2).\n');

f = zeros(6,1);
B = zeros(6,2);

for i=1:3,
  k1 = 2*i - 1;
  k2 = k1 + 1;
  
  f(k1:k2,1) = -x(:,i);
  B(k1:k2,:) = [-y(1,i), -y(2,i); -y(2,i), y(1,i)];
  
end

N = B'*W*B;
iN = inv(N);
t = B'*W*f;

delta = iN * t;

c = delta(1);
d = delta(2);

fprintf('Then c = %.6f, d = %.6f\n', c, d)
fprintf('a = %.6f, b = %.6f\n', c/(c^2+d^2), d/(c^2+d^2))