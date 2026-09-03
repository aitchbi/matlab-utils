function t = hb_runcmd(cmd,errmsg)
tic;
[sts,log] = system(cmd);
t = toc;
if sts==0
    return;
end
sprintf('*** system run command log: \n\n');
log %#ok<NOPRT>
fprintf('\n.%s [run time: %d minutes]', errmsg, round(t/60));
error(errmsg)
end
