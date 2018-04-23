function [mu, P, y_prefit, y_postfit] = KF(mu, P, F, H, Q, R, z)
    % predict
    mu_predict = F*mu;
    P_predict = F*P*F' + Q;
    % update
    y_prefit = z - H*mu_predict;
    K = P_predict*H'/(R+H*P_predict*H');
    mu = mu_predict + K*y_prefit;
    P = (eye(length(mu)) - K*H)*P_predict*(eye(length(mu)) - K*H)' + ...
        K*R*K';
    y_postfit = z - H*mu;
