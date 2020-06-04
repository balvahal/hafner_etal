function matrix = smoothMatrix(matrix, span, method)
    for i=1:size(matrix,1)
        if(strcmp(method, 'median'))
            trace = medfilt1(matrix(i,matrix(i,:) ~= -1), span);
        else
            trace = smoothdata(matrix(i,matrix(i,:) ~= -1), method, span);
        end
        matrix(i,matrix(i,:) ~= -1) = trace;
    end
end