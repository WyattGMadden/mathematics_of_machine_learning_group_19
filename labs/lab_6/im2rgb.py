%% load and convert an IMage file 2 a sequence of RGB triples
%%  Load an image as a sequence of RGB triples.
%% Inputs:
%%  ImageFile - A string with the path to (name of) an image file to load.
%% Outputs:
%%  X - 3-by-(dims(0)*dims(1)) matrix with an RGB triple in each column.
%%  I - dims(0)-by-dims(1)-by-dims(2) matrix of integers with the raw image data.
%%  dims - A triple of numbers with the height, width, and number of channels (3 for RGB).

function [X,I,dims] = im2rgb(ImageFile)
    I = imread(ImageFile);

    W = size(I,2);
    H = size(I,1);
    XI = cast(reshape(I,[W*H,3]),'double')/255.0;
    dims = size(I);
    X = XI';
end