clc;clear;
[filepath]=fileparts([mfilename('fullpath'),'.m']); 
addpath(genpath(filepath))
cd(filepath)

results = [];
count = 0;

tol = 1e-8;

for x1 = 1:49
    for x2 = 1:49
        val1 = (3.8303 * x1 - 8.8764 * x2) / (3.8303 * x1);
        if val1 <= 0.9 + tol || val1 >= 1.1 - tol
            continue;
        end
        for x3 = 1:49
            for x4 = 1:49
                val2 = (3.8303 * x3 - 8.6925 * x4) / (3.8303 * x3);
                if val2 > 0.9 + tol && val2 < 1.1 - tol
                    results(end+1, :) = [x1, x2, x3, x4];
                    count = count + 1;
                    if count >= 10
                        break;
                    end
                end
            end
            if count >= 10
                break;
            end
        end
        if count >= 10
            break;
        end
    end
    if count >= 10
        break;
    end
end

disp('前10组更严格匹配 Mathematica 条件的整数解：');
disp(results);



