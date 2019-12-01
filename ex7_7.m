%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 7.7 which solves 
% Example 6.1 using the technique of indirect observations.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 163-164.
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

% Table 6.1: observation value and standard deviation 
obs = [ ...
16.5, 0.1;   % 1 (mm)
 3.8, 0.1;   % 2 (mm)
20.4, 0.1;   % 3 (mm)
10.0, 0.05;  % 4 (m)
 8.0, 0.05]; % 5 (m) 

% principal distance
p = 100; % mm

% parameter vector definitionin meters:
% x(1) = x_1 coordinate of P 
% x(2) = x_2 coordinate of P 
% x(3) = distance between S_1 and S_2
% x(4) = distance between S_2 and S_3

% nonliner funional model: 
% l_1 + v_1 = p * x_1 / x_2
% l_2 + v_2 = p * (x_3 - x_1)/x_2
% l_3 + v_3 = p * (x_3 + x_4 - x_1)/x_2
% l_4 + v_4 = x_3
% l_5 + v_5 = x_4
% 
% From 
%  l_1 + l_2 + v_1 + v_2 = p * x_3 / x_2
% we can initialie
%  x0(2) = p * l_4 / (l_1 + l_2)
%  x0(1) = l_1 * x0(2) / p = l_1 * l_4 /(l_1 + l_2)

% initial value of the parameters
x0 = [ ...
obs(1,1)*obs(4,1)/(obs(1,1)+obs(2,1));
p*obs(4,1)/(obs(1,1)+obs(2,1));
obs(4,1);
obs(5,1)];

% linearzed model
B = [...
-p/x0(2),               p*x0(1)/x0(2)^2,        0,        0;
 p/x0(2),       p*(x0(3)-x0(1))/x0(2)^2, -p/x0(2),        0;
 p/x0(2), p*(x0(3)+x0(4)-x0(1))/x0(2)^2, -p/x0(2), -p/x0(2); 
       0,                             0,       -1,        0;
	   0,                             0,        0,       -1];
% with reference variance = 0.01
ref_var = 0.01;
W = diag([1,1,1,4,4]);

% normal matrix
N = B' * W * B

% 
f = zeros(5,1);
f(1) = p*x0(1)/x0(2)-obs(1,1);
f(2) = p*(x0(3)-x0(1))/x0(2)-obs(2,1);
f(3) = p*(x0(3)+x0(4)-x0(1))/x0(2)-obs(3,1);
f(4) = x0(3)-obs(4,1);
f(5) = x0(4)-obs(5,1);

t = B' * W * f

Qx = inv(N)

delta = Qx * t

% Adjusted parameters
x = x0 + delta
