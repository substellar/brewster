# settings.py

# def init():
#     global runargs
#     global linelist
#     global samples
#     samples = []
#     runargs = 0


def init(args_instance=None):
    global runargs
    global linelist
    global samples
    # global saverunargs
    samples = []
    if args_instance is not None:
        runargs = args_instance
    else:
        runargs = 0

    

        # saverunargs = (
        #     args_instance.gases_myP,
        #     args_instance.chemeq,
        #     args_instance.dist,
        #     args_instance.cloudtype,
        #     args_instance.do_clouds,
        #     args_instance.gaslist,
        #     args_instance.gasnames,
        #     args_instance.gasmass,        
        #     args_instance.cloudnum,
        #     args_instance.inlinetemps,
        #     args_instance.coarsePress,
        #     args_instance.press,
        #     args_instance.inwavenum,
        #     args_instance.linelist,
        #     args_instance.cia,
        #     args_instance.ciatemps,
        #     args_instance.use_disort,
        #     args_instance.fwhm,
        #     args_instance.obspec,
        #     args_instance.proftype,
        #     args_instance.do_fudge,
        #     args_instance.prof,
        #     args_instance.do_bff,
        #     args_instance.bff_raw,
        #     args_instance.ceTgrid,
        #     args_instance.metscale,
        #     args_instance.coscale
        # )
