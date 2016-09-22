pro read_para,dir, para

    p = Python.Run("import HybridParams")
    p = Python.Run("p = HybridParams.HybridParams("+string(34B)+dir+string(34B)+")")
    p = Python.p
    p = p.para
    para = p.ToStruct()
    
return
end
