%% RGB colors 2 IMage conversion
%%  Converts a sequence of RGB colors representing an image with the given
%%   dimensions into a MATLAB image.
%% Inputs:
%%  X - A 3-by-(dims(0)*dims(1)) matrix where each column represents the RGB
%%    values of a pixel in the image.
%%  dims - Three integers representing the size of the image to create.
%% Outputs:
%%  I - A dims(0)-by-dims(1)-by-3 matrix of integers representing the RGB image.

function [I] = rgb2im(X,dims)
	if (numel(dims) ~= 3), error('rgb2im requires dims to be three integer values.'); end
	if (dims(1)*dims(2) ~= size(X,2)), error('rgb2im number of RGB triples does not match image size.'); end
	if (size(X,1) ~= 3), error('rgb2im must have X with three rows (RGB values).'); end
    I = cast(255.0*reshape(X',[dims(1),dims(2),3]),'uint8');
end