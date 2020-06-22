function [q_new,qdot_new] = finitediff(w_dot,q_dot,q)
delta_t=0.01;

qdot_new = q_dot + w_dot*delta_t;
q_new =    q + qdot_new *delta_t;

end