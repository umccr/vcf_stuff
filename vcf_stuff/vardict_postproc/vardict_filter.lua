function check_msi(msi, msilen, af)
    if (msi == nil or msilen == nil or af == nil) then
        return false
    end
    local change_len = math.abs(ref.len() - alt.len())
    local msi_fail = false
    if (change_len == 1 and msi > 1) then
        if (msi <= 2 and af < 0.005) then
            msi_fail = true
        end
    end
    return msi_fail
end
