%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 9.1 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 219-221
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

% Known coordinates: P(0,5), O(0,0), Q(5,0)

% distance observations: RP, RO, RQ
L = [ 5.10; 7.07; 4.90];

% approximate values for the unknown coordinates of R

x1 = 5; x2 = 5;

% Linearize the distance measurements to the form: v + B delta = f

L0 = [ sqrt(x1^2+(x2-5)^2); sqrt(x1^2 + x2^2); sqrt((x1-5)^2+x2^2)];

B = [...
-x1/L0(1), -(x2-5)/L0(1);
-x1/L0(2), -x2/L0(2);
-(x1-5)/L0(3), -x2/L0(3)];

f = L0 - L;

N = B' * B;
t = B' * f;

invN = inv(N);
delta = invN * t;

% Given the constraint x2 = 5, compute the corresponding value of x1
% C delta = g, where C = [0 1]

C = [0 1];
g = 0;

M = C * invN * C';

kc = inv(M)*(g - C * delta);

dd = invN * C' * kc;

% Compute delta with the constraint
delta = delta + dd;

