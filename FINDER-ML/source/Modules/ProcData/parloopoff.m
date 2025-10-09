function parloopoff(parameters)

if parameters.parallel.on == 1
    delete(gcp('nocreate'));
end

end
