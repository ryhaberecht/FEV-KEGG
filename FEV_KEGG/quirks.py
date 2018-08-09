"""
Sometimes, things do not work as expected. This is when quirk workarounds save the day.

Warnings
--------
These quirks are highly time-specific, they can change any moment! Be sure to re-evaluate the workarounds regularly.
"""

# some organisms are listed, but do not return anything useful when queried via the REST API
NON_EXISTING_ORGANISMS = ['lni', 'scla', 'pavl', 'our', 'rox', 'dei', 'miq', 'phz', 'mela', 'simp', 'grs', 'deo', 'mdv', 'hag', 'ema', 'pset', 'stro', 'bmur', 'malk', 'lyb', 'ptc', 'pamg', 'vit', 'syo', 'thas', 'paih', 'otk', 'laca', 'pio', 'fsa', 'lcy', 'mgg', 'psai', 'bzg', 'cbae', 'git', 'aue', 'mbas', 'kit', 'vdb', 'aid', 'smur', 'phr', 'melm', 'mee'] # as of: 2018-04-06
"""
Some organisms in KEGG supposedly exist, but when retrieved, their data returns empty, which raises an Error.

Warnings
--------
Today, some organisms might be in this list without reason, others might be missing.
Please try to fetch the organisms in this list by executing

::

    org = FEV_KEGG.Organism.Organism('lni')
    org.getPathways()

If this fails, then lni still belongs in this list.
"""

# there is no naming convention by which to differentiate metabolic pathways from structural, signalling, etc. pathways
METABOLIC_PATHWAYS = ['01100', '01110', '01120', '01130', '01200', '01210', '01212', '01230', '01220', '00010', '00020', '00030', '00040', '00051', '00052', '00053', '00500', '00520', '00620', '00630', '00640', '00650', '00660', '00562', '00190', '00195', '00196', '00710', '00720', '00680', '00910', '00920', '00061', '00062', '00071', '00072', '00073', '00100', '00120', '00121', '00140', '00561', '00564', '00565', '00600', '00590', '00591', '00592', '01040', '00230', '00240', '00250', '00260', '00270', '00280', '00290', '00300', '00310', '00220', '00330', '00340', '00350', '00360', '00380', '00400', '00410', '00430', '00440', '00450', '00460', '00471', '00472', '00473', '00480', '00510', '00513', '00512', '00515', '00514', '00532', '00534', '00533', '00531', '00563', '00601', '00603', '00604', '00540', '00550', '00511', '00730', '00740', '00750', '00760', '00770', '00780', '00785', '00790', '00670', '00830', '00860', '00130', '00900', '00902', '00909', '00904', '00906', '00905', '00981', '00908', '00903', '00281', '01052', '00522', '01051', '01059', '01056', '01057', '00253', '00523', '01054', '01053', '01055', '00940', '00945', '00941', '00944', '00942', '00943', '00901', '00403', '00950', '00960', '01058', '00232', '00965', '00966', '00402', '00311', '00332', '00261', '00331', '00521', '00524', '00525', '00231', '00401', '00404', '00405', '00333', '00254', '00362', '00627', '00364', '00625', '00361', '00623', '00622', '00633', '00642', '00643', '00791', '00930', '00351', '00363', '00621', '00626', '00624', '00365', '00984', '00980', '00982', '00983', '00970']
"""
KEGG PATHWAY does not provide a reliable naming convention for pathways which are part of metabolism and actually do contain 'Compounds' and Enzymes/Genes.

This is the list of all pathway names which currently provide usable information about metabolism.

Warnings
--------
It may change any time! Compare with `KEGG PATHWAY <http://www.kegg.jp/kegg/pathway.html>`_ to see whether you want to categorise different pathways as 'metabolic'.
"""