function IR = rotationMatrix(phi, d)

    theta = phi*ones(d-1,1);
    % theta = [zeros(d-2,1); phi]; % hier wird nur der letzte Winkel gedreht
    I=eye(d);
    IR=I;
    vv=I(d,:);
    for r=d-1:-1:1
        uu=I(r,:);
        IR=IR*RotMatrix(-theta(r),uu,vv);
    end

end