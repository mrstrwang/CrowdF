from random import Random

from operations.models import EmailVerifyCode
from django.core.mail import send_mail, EmailMessage
from CrowdFunding.settings import EMAIL_FROM


def random_str(random_length=8):
    str = ''
    # 生成字符串的可选字符串
    chars = 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz0123456789'
    length = len(chars) - 1
    random = Random()
    for i in range(random_length):
        str += chars[random.randint(0, length)]
    return str


def send_email_code(email, send_type="register", user_id=None):
    resend = EmailVerifyCode.objects.filter(email=email, send_type=send_type)
    if resend:
        resend.delete()
    email_record = EmailVerifyCode()
    if send_type == "verify_id":
        code = random_str(6)
    else:
        code = random_str(16)
    email_record.verify_code = code
    email_record.email = email
    email_record.send_type = send_type
    email_record.save()

    email_title = ""
    email_body = ""

    if send_type == "register":
        email_title = "燎原众筹网 注册激活链接"
        email_body = "欢迎注册燎原众筹网:  请点击下面的链接激活你的账号: http://192.168.20.37:8000/users/active/{0}LyZcW{1}".format(code, user_id)
        msg = EmailMessage(email_title, email_body, EMAIL_FROM, [email])
        msg.content_subtype = "html"
        send_status = msg.send()
        if send_status:
            pass

    elif send_type == "verify_id":
        email_title = "燎原众筹网 实名认证验证码"
        email_body = "燎原众筹网:  你的邮箱验证码为: {0}".format(code)
        msg = EmailMessage(email_title, email_body, EMAIL_FROM, [email])
        msg.content_subtype = "html"
        send_status = msg.send()
        if send_status:
            pass
