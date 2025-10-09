function X = wipeTrainingSub(X)

if ~isstruct(X)
    return

elseif isstruct(X)

    F = fields(X);
    for i = 1:length(F)
        
        field = F{i};

        if contains(field, 'train', 'IgnoreCase', true)
            X = rmfield(X, field);
        else
            X.(field) = wipeTrainingSub(X.(field));
        end


    end

end





end