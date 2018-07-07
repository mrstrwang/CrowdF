from django.urls import path
from users.views import *

app_name = 'users'

urlpatterns = [
    path('login/', LoginView.as_view(), name='login'),
    path('logout/', LogoutView.as_view(), name='logout'),
    path('register/', RegisterView.as_view(), name='register'),
    path('active/<str:active_code>/', ActiveView.as_view(), name='active'),
    path('user_center/', UserCenterView.as_view(), name='user_center'),
    path('user_crowdfunding/', UserCrowdFundingView.as_view(), name='user_crowdfunding'),
    path('user_accttype/', UserAcctTypeView.as_view(), name='user_accttype'),
]