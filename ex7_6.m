%% Eun Hwan Shin, 2019
%
% This file is the Octave implementation of Example 7.6 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 161-162
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

% observed angles (deg,min,sec) of a simple plane triangle.
obs = [ ...
40, 19, 02;
70, 30, 01;
69, 11, 00];

% angles in arc seconds
o_deg = obs(:,1) + obs(:,2)/60 + obs(:,3)/3600;

B = [ ...
-1,  0;
 0, -1;
 1,  1];

f = [ ...
-o_deg(1);
-o_deg(2);
180 - o_deg(3)];

N = B' * B;
t = B' * f;

Qxx = inv(N);
delta = Qxx * t;

x_dms = deg2dms(delta)

v = f - B * delta;

fprintf('Residuals in arc-sec\n')
v_as = v*3600

%Cofactor matrix of the adjusted observation
Qll = B * Qxx * B'