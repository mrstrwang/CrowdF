from django.urls import path

from orders.views import *

app_name = 'orders'

urlpatterns = [
    path('pay_step/<int:item_id>/', PayStep1View.as_view(), name='pay_step1'),
    path('pay_step2/', PayStep2View.as_view(), name='pay_step2'),
    path('pay_order/', PayOrderView.as_view(), name='pay_order'),
    path('pay_money/', PayMoneyView.as_view(), name='pay_money'),
    path('pay_success/', PaySuccessView.as_view(), name='pay_success'),
]
