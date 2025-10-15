function f_func = hb_fsl_get_func(d_fsl, func)
f_func = fullfile(d_fsl, 'bin', func);
d = exist(f_func, 'file');
switch d
    case 2
        % all good
    case 7
        f_func = fullfile(d_fsl, 'bet'); % maybe f_fsl is the bin path?
        assert(exist(f_func, 'file')==2);
    otherwise
        error('fishy FSL path..');
end
end
