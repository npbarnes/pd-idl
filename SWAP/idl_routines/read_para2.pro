pro read_para,dir, para

    p = Python.Run("import HybridParams")
    p = Python.Run("p = HybridParams.HybridParams("+string(34B)+dir+string(34B)+")")
    p = Python.p
    para = p.para
    
return
end
