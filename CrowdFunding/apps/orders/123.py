from alipay import AliPay
from CrowdFunding.settings import ALIPAY_URL, BASE_DIR

def get_alipay():
    app_private_key_string = open(BASE_DIR + "/utils/app_private_key.pem").read()
    alipay_public_key_string = open(BASE_DIR + "/utils/app_public_key.pem").read()
    alipay = AliPay(
        appid="2016091500515408",
        app_notify_url='http://127.0.0.1:8000/orders/pay_success/?order=9',  # 默认回调url
        app_private_key_string=app_private_key_string,
        alipay_public_key_string=alipay_public_key_string,  # 支付宝的公钥，验证支付宝回传消息使用，不是你自己的公钥,
        sign_type="RSA2",  # RSA 或者 RSA2
        debug=True  # 默认False
    )
    return alipay


def get_order_string():
    subject = "测试订单"
    order_string = get_alipay().api_alipay_trade_page_pay(
        out_trade_no="20161112",
        total_amount=0.01,
        subject=subject,
        return_url='http://127.0.0.1:8000/orders/pay_success/?order=9',
        notify_url='http://127.0.0.1:8000/orders/pay_success/?order=8'   # 可选, 不填则使用默认notify url
    )
    pay_url = ALIPAY_URL + '?' + order_string
    print(pay_url)
    return pay_url


if __name__ == '__main__':
    alipay = get_alipay()
    result = alipay.api_alipay_trade_query(trade_no='2018051221001004950200286692')
    print(result)