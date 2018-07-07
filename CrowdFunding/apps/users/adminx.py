import xadmin
from xadmin import views

from users.models import UserProfile, UserAddress


class BaseSetting(object):
    enable_themes = True
    use_bootswatch = True


class GlobalSettings(object):
    site_title = "燎原众筹平台"
    site_footer = "燎原众筹平台"


xadmin.site.unregister(UserProfile)
xadmin.site.register(UserProfile)
xadmin.site.register(UserAddress)


xadmin.site.register(views.BaseAdminView, BaseSetting)
xadmin.site.register(views.CommAdminView, GlobalSettings)