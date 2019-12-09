%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 7.11.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 170-172
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

table7_1

% Weight matrix
W = diag((0.1*dist).^-1);

% v + B x = f, where x consts of the elevations of points B, C, D, and E.
B = [ ...
-1,  0,  0,  0;
 1, -1,  0,  0;
 0,  1,  0,  0;
 1,  0, -1,  0;
 0,  0,  1, -1;
 0, -1,  0,  1;
 0,  0,  0,  1;
 0,  1, -1,  0];
 
 f = -dE;
 
 N = B' * W * B;
 t = B' * W *f;
 invN = inv(N);
 
 delta = invN * t;
 
 x = delta + 800;
 
 v = f - B * delta;
 
 ref_var =  v'*W*v / 4;
 
 % covariance of the elevations
 Cxx = ref_var * invN;