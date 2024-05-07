function hb_runcmd(cmd,errmsg)
[sts,log] = system(cmd);
if sts==0
    return;
end
sprintf('*** system run command log: \n\n');
log %#ok<NOPRT>
error(errmsg)
end
