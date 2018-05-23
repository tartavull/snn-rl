# from http://stackoverflow.com/questions/18019477/how-can-i-play-a-local-video-in-my-ipython-notebook
import io;import base64;from IPython.display import HTML

def playVideo(videoName):
	video = io.open(videoName, 'r+b').read();encoded = base64.b64encode(video)
	return HTML(data='''<video alt="test" controls><source src="data:video/mp4;base64,{0}" type="video/mp4" /></video>'''.format(encoded.decode('ascii')))