function isapprox_act_profs(act_prof1, act_prof2)
    length(act_prof1) == length(act_prof2) || return false
    for (action1, action2) in zip(act_prof1, act_prof2)
        isapprox(action1, action2) || return false
    end
    return true
end

function isapprox_vecs_act_profs(vec_act_profs1, vec_act_profs2)
    length(vec_act_profs1) == length(vec_act_profs2) || return false
    for act_prof in vec_act_profs1
        any(x -> isapprox_act_profs(x, act_prof),
            vec_act_profs1) || return false
    end
    return true
end
