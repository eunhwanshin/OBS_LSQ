%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 4.3 in
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 87.
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

% Given [y_1; y_2] = A [x_1; x_2], z = x_1 + x_2, where
A = [2, 1; 1, 2];
Cyy = [3, 1; 1, 3];

% Compute sigma_z^2, Czy, and Cxy

Jxy = inv(A);

Cxx = Jxy * Cyy * Jxy';

Jzx = [1, 1];
Czz = Jzx * Cxx * Jzx';

fprintf('sigma_z^2 = %.6f\n', Czz);

Czy = Jzx * Jxy * Cyy;

Cxy = Jxy * Cyy;

% Standard ellipsoid for x
[v, lambda] = eig(Cxx);
ang = atan(v(2,2)/v(1,2))*180/pi;
fprintf('Standard ellipsoid for x\n');
fprintf('semi-major = %.6f, semi-minor = %.6f\n', sqrt(lambda(2,2)), sqrt(lambda(1,1)))
fprintf('orientation of semi-major = %.6f\n', ang);

% Standard ellipsoid for z: one dimentional case
std_z = sqrt(Czz);
fprintf('Standard ellipsoid for z = [-%.6f, %.6f]\n', std_z, std_z);



