%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 9.3, which solves
% the problem in Example 9.2 by elimination of the constraint. 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 231-232.
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

% condition equation: v + B delta = f
B = [ ...
 -3,  1, -1;
  1, -4, -2];
  
f = [2; 18]; 

% contraint equation: C delta = g
C = [...
 -1,  1, 1;
  1, -3, 1];
g = [ 2; 26];

% a priori weight matrix
W = diag([2, 0.5]);

% Direct solution: since the number of condition equations is smaller than
% that of the unknowns N is singular.

N = B' * W * B;
t = B' * W * f;

% However, the augmented matrix K is not singular.
K = [-N, C'; C, zeros(2)];

invK = inv(K);

x = invK * [-t; g];
delta = x(1:3)
Qdd = -invK(1:3,1:3);

% Solution by Elimination:

% partition delta = [delta1; delta2] and thus
B1 = B(:,1:2);
B2 = B(:,3);

C11 = C(:,1:2);
C12 = C(:,3);

invC11 = inv(C11)

B_bar = -B1 * invC11 * C12 + B2;
f_bar = f - B1 * invC11 * g;

N_bar = B_bar' * W * B_bar;
t_bar = B_bar' * W * f_bar;

delta2 = t_bar/N_bar;
delta1 = invC11 * (g - C12 * delta2);

fprintf('delta - [delta1; delta2] = \n')
delta - [delta1; delta2]
