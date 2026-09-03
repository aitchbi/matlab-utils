function Y = hb_rs_timecourse_bpfilter(X, band, TR)
% HB_RS_TIMECOURSE_BPFILTER  band-pass filter fMRI time courses.
%
% this function is based on the temporal filtering implementation in
% conn_filter.m from the CONN functional connectivity toolbox [1].
%
% inputs
%   X   : T x P matrix of time courses
%           T = number of time points
%           P = number of parcels/voxels
%
%   band: [low high] frequency band in Hz, e.g. [0.01 0.1]
%
%   TR  : repetition time in seconds
%
% output
%   Y   : T x P band-pass filtered time courses
%
% [1] Whitfield-Gabrieli S, Nieto-Castanon A. Conn: a functional
% connectivity toolbox for correlated and anticorrelated brain networks.
% Brain Connectivity. 2012;2(3):125-141.

    T = size(X,1);

    fy = fft([X; flipud(X)], [], 1);

    f = 0:size(fy,1)-1;
    f = min(f, size(fy,1)-f);

    idx = f < band(1)*(TR*size(fy,1)) | ...
          f >= band(2)*(TR*size(fy,1));

    fy(idx,:) = 0;

    Y = real(ifft(fy, [], 1));
    Y = Y(1:T,:);

    assert_bpfiltered(Y, band, TR);
end

%==========================================================================
function assert_bpfiltered(y, band, TR)
% verify that a CONN-style band-pass filtered time course contains no
% appreciable frequency components outside the requested frequency band.

    fy = fft([y; flipud(y)], [], 1);

    f = 0:size(fy,1)-1;
    f = min(f, size(fy,1)-f);

    idx_out = f < band(1)*(TR*size(fy,1)) | ...
              f >= band(2)*(TR*size(fy,1));

    % relative magnitude of largest out-of-band Fourier coefficient
    max_all = max(abs(fy), [], 1);
    max_out = max(abs(fy(idx_out,:)), [], 1);

    rel_out = max_out ./ max(max_all, eps);
    
    msg = 'out-of-band signal remains.';
    
    assert(all(rel_out < 1e-10), ...
        sprintf('fishy: band-pass filtering verification failed: %s', msg));
end