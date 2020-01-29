function A = LowerFrameBound(X)
    d = size(X,1);
    S = zeros(d,d);
    for i=1:d
        Sd = sum(X*diag(X(i,:)),2);
        S(i,:) = Sd;
    end

    A = min(eig(S));
end
