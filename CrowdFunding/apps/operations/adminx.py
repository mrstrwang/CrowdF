import xadmin

from operations.models import UserSupport, UserInterest, Banner, Advertise, EmailVerifyCode


xadmin.site.register(UserSupport)
xadmin.site.register(UserInterest)
xadmin.site.register(Banner)
xadmin.site.register(Advertise)
xadmin.site.register(EmailVerifyCode)
