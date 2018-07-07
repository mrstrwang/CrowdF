import time
from datetime import timedelta, datetime

from django.shortcuts import render
from django.views.generic.base import View
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect

from projects.models import ProjectItem, Project
from users.models import UserAddress
from orders.models import UserPayOrder
from alipay import AliPay
from CrowdFunding.settings import ALIPAY_URL, BASE_DIR

# Create your views here.


class PayStep1View(View):
    def get(self, request, item_id):
        item = ProjectItem.objects.get(id=item_id)
        item.deadline_time = item.project.create_time + timedelta(days=item.project.deadline)
        return render(request, 'pay-step-1.html', {'item': item, })


class PayStep2View(View):
    def get(self, request, item_id):
        item = ProjectItem.objects.get(id=item_id)
        item.deadline_time = item.project.create_time + timedelta(days=item.project.deadline)
        return render(request, 'pay-step-2.html', {'item': item, })

    def post(self, request):
        item_id = request.POST.get('item_id', '')
        nums = request.POST.get('nums', '')
        user = request.user
        user_address_list = UserAddress.objects.filter(user_id=user.id)
        if request.POST.get('new_address', ''):
            user_address = UserAddress.objects.filter(user_id=user.id)
            # name = request.POST.get('name', '')
            # mobile = request.POST.get('mobile', '')
            # address = request.POST.get('address', '')
            name, mobile, address = [request.POST.get(x) for x in ['name', 'mobile', 'address']]
            if not user_address:
                UserAddress.objects.create(name=name, mobile=mobile, address=address, user_id=user.id, is_default=True)
            elif len(user_address) > 4:
                return HttpResponse('{"status":"fail"}', content_type='application/json')
            else:
                UserAddress.objects.create(name=name, mobile=mobile, address=address, user_id=user.id)
                return HttpResponse('{"status":"success"}', content_type='application/json')

        item = ProjectItem.objects.get(id=item_id)
        item.total_pay = '%.2f' % (float(item.item_value) * float(nums))

        return render(request, 'pay-step-2.html', {'user_address_list': user_address_list,
                                                   'item': item,
                                                   'nums': nums, })


class PayOrderView(View):
    def get(self, request):
        return render(request, 'index.html')

    def post(self, request):
        item_id = request.POST['item_id']
        item = ProjectItem.objects.get(id=item_id)
        order = UserPayOrder()
        order.item_nums = request.POST['order_nums']
        order.freight = item.item_freight
        order.total_pay = request.POST['total_pay']
        order.need_invoice = request.POST['order_need_invoice']
        order.info_invoice = request.POST['order_invoice_info']
        order.remarks = request.POST['order_remarks']
        order.address_id = request.POST['order_address']
        order.item_id = request.POST['item_id']
        order.user_id = request.user.id
        order.order_status = 'waiting'
        order.save()
        order.order_no = '2018' + str(int(time.time())) + str(item_id) + 'LyZc' + str(order.id)
        order.save()
        return JsonResponse({'status': 'success', 'order_id': order.id})


class PayMoneyView(View):
    def get(self, request):
        order_id = int(request.GET.get('order'))
        order = UserPayOrder.objects.filter(id=order_id)
        if order:
            order = order[0]
            subject = order.item.project.name + order.item.item_name
            out_trade_no = order.order_no
            total_amount = str(order.total_pay)
            pay_url = get_order_string(subject, out_trade_no, total_amount)
            return HttpResponseRedirect(pay_url)
        return HttpResponse('订单信息错误！', content_type='text/html')


class PaySuccessView(View):
    def get(self, request):
        alipay = get_alipay()
        trade_no = request.GET.get('trade_no')
        result = alipay.api_alipay_trade_query(trade_no=trade_no)
        code = result.get('code')
        print(result)
        while True:
            if code == '10000' and result.get('trade_status') == 'TRADE_SUCCESS':
                order_id = request.GET.get('out_trade_no', '').split('LyZc')[1]
                order = UserPayOrder.objects.get(id=order_id)
                order.order_status = 'paid'
                order.pay_order_no = request.GET.get('trade_no')
                order.pay_time = datetime.now()
                order.save()
                project_id = order.item.project.id
                project = Project.objects.get(id=project_id)
                project.support_fund += order.total_pay
                project.support_nums += 1
                project.save()
                return render(request, 'pay_success.html')
            elif code == '40004' or (code == '10000' and result.get('trade_status') == 'WAIT_BUYER_PAY'):
                time.sleep(5)
                continue
            else:
                return HttpResponse('出错了')




def get_alipay():
    app_private_key_string = open(BASE_DIR + "/utils/app_private_key.pem").read()
    alipay_public_key_string = open(BASE_DIR + "/utils/app_public_key.pem").read()
    alipay = AliPay(
        appid="2016091500515408",
        app_notify_url=None,  # 默认回调url
        app_private_key_string=app_private_key_string,
        alipay_public_key_string=alipay_public_key_string,  # 支付宝的公钥，验证支付宝回传消息使用，不是你自己的公钥,
        sign_type="RSA2",  # RSA 或者 RSA2
        debug=True  # 默认False
    )
    return alipay


def get_order_string(subject, out_trade_no, total_amount):
    subject = subject
    order_string = get_alipay().api_alipay_trade_page_pay(
        out_trade_no=out_trade_no,
        total_amount=total_amount,
        subject=subject,
        return_url='http://127.0.0.1:8000/orders/pay_success/',
    )
    pay_url = ALIPAY_URL + '?' + order_string
    print(pay_url)
    return pay_url
